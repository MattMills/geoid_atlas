//! Recursive Gaussian state estimation and fixed-interval historical smoothing.
//!
//! State order is [x, y, z, vx, vy, vz] in one declared frame/time coordinate.
//! Frames are deterministic; clock, ephemeris, and sensor systematic uncertainty
//! are NOT estimated by this six-state model. Measurements must be conditionally
//! independent or decorrelated by a caller-supplied measurement model.
//!
//! ```
//! use geoid_atlas::estimation::*;
//! use geoid_atlas::frames::{Epoch, State, Vec3};
//! let initial = TrackState::new(Epoch::new(0.0)?,
//!     State::stationary(Vec3::ZERO), Covariance6::isotropic(1.0, 1.0)?)?;
//! let mut tracker = RecursiveTracker::new(initial, Vec3::ZERO)?;
//! let observation = LinearMeasurement {
//!     coefficients: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
//!     observed: 1.0, noise_std: 1.0,
//! };
//! let record = tracker.assimilate("sensor-1/sample-1", Epoch::new(1.0)?,
//!                                 &observation, Some(9.0))?;
//! assert!(record.diagnostics.accepted);
//! let history = smooth_history(&[initial, record.filtered], &[record.predicted])?;
//! assert!((history[0].mean.position_m.x - 1.0 / 3.0).abs() < 1e-12);
//! # Ok::<(), geoid_atlas::Error>(())
//! ```
use crate::frames::{Epoch, FrameGraph, FrameId, Pose, State, Vec3};
use crate::spacetime::{ConstantAcceleration, HermiteWorldline, Worldline, bistatic_path};
use crate::{Error, Result, finite};
use std::collections::HashSet;

pub type Matrix6 = [[f64; 6]; 6];
pub type Vector6 = [f64; 6];

fn invalid(message: &str) -> Error {
    Error::InvalidInput(message.into())
}
fn identity() -> Matrix6 {
    let mut result = [[0.0; 6]; 6];
    for (i, row) in result.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    result
}
fn transpose(a: Matrix6) -> Matrix6 {
    let mut result = [[0.0; 6]; 6];
    for (i, row) in result.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = a[j][i];
        }
    }
    result
}
fn multiply(a: Matrix6, b: Matrix6) -> Matrix6 {
    let mut result = [[0.0; 6]; 6];
    for (i, row) in result.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = (0..6).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    result
}
fn add(a: Matrix6, b: Matrix6) -> Matrix6 {
    let mut result = a;
    for (i, row) in result.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value += b[i][j];
        }
    }
    result
}
fn subtract(a: Matrix6, b: Matrix6) -> Matrix6 {
    let mut result = a;
    for (i, row) in result.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value -= b[i][j];
        }
    }
    result
}
fn multiply_vector(a: Matrix6, x: Vector6) -> Vector6 {
    std::array::from_fn(|i| (0..6).map(|j| a[i][j] * x[j]).sum())
}
fn state_vector(state: State) -> Vector6 {
    [
        state.position_m.x,
        state.position_m.y,
        state.position_m.z,
        state.velocity_m_s.x,
        state.velocity_m_s.y,
        state.velocity_m_s.z,
    ]
}
fn vector_state(x: Vector6) -> Result<State> {
    Ok(State {
        position_m: Vec3::new(x[0], x[1], x[2])?,
        velocity_m_s: Vec3::new(x[3], x[4], x[5])?,
    })
}
fn transition(dt: f64) -> Matrix6 {
    let mut f = identity();
    for (i, row) in f.iter_mut().take(3).enumerate() {
        row[i + 3] = dt;
    }
    f
}

// Operate on correlations for validation/solving to handle mixed position and
// velocity units without treating small velocity variances as numerically zero.
fn correlation(matrix: Matrix6) -> Result<(Matrix6, Vector6)> {
    let mut scales = [0.0; 6];
    for (i, value) in scales.iter_mut().enumerate() {
        finite(matrix[i][i], "covariance diagonal")?;
        if matrix[i][i] < 0.0 {
            return Err(invalid("covariance diagonal must be nonnegative"));
        }
        *value = matrix[i][i].sqrt();
    }
    let mut normalized = [[0.0; 6]; 6];
    for (i, row) in normalized.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            finite(matrix[i][j], "covariance entry")?;
            if scales[i] == 0.0 || scales[j] == 0.0 {
                if matrix[i][j] != 0.0 {
                    return Err(invalid("zero-variance covariance row must be zero"));
                }
            } else {
                *value = matrix[i][j] / scales[i] / scales[j];
                finite(*value, "covariance correlation")?;
            }
        }
    }
    Ok((normalized, scales))
}
fn cholesky(matrix: Matrix6, allow_semidefinite: bool) -> Result<Matrix6> {
    let mut l = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let residual = matrix[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            finite(residual, "covariance factorization")?;
            if i == j {
                if residual < -1e-12 || (!allow_semidefinite && residual <= 1e-14) {
                    return Err(invalid(
                        "covariance must be positive semidefinite; solve requires positive definite covariance",
                    ));
                }
                l[i][j] = residual.max(0.0).sqrt();
            } else if l[j][j] == 0.0 {
                if residual.abs() > 1e-12 {
                    return Err(invalid("covariance is not positive semidefinite"));
                }
            } else {
                l[i][j] = residual / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Symmetric positive-semidefinite full Cartesian state covariance.
/// Diagonal units are m² for positions and (m/s)² for velocities;
/// off-diagonal units are the corresponding products. Exact singular covariances
/// are allowed for filtering, but historical smoothing requires invertible predictions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Covariance6 {
    matrix: Matrix6,
}
impl Covariance6 {
    pub fn new(matrix: Matrix6) -> Result<Self> {
        let (normalized, _) = correlation(matrix)?;
        for (i, row) in normalized.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if (*value - normalized[j][i]).abs() > 1e-12 {
                    return Err(invalid("covariance must be symmetric"));
                }
            }
        }
        let mut symmetric = matrix;
        for (i, row) in symmetric.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = 0.5 * matrix[i][j] + 0.5 * matrix[j][i];
            }
        }
        cholesky(correlation(symmetric)?.0, true)?;
        Ok(Self { matrix: symmetric })
    }
    // Products such as J P J^T are symmetric analytically. Remove floating-point
    // antisymmetry BEFORE validating their normalized covariance, particularly
    // when smoothing subtracts large predicted terms. User inputs still go
    // through new() and are checked for symmetry without this correction.
    fn calculated(matrix: Matrix6) -> Result<Self> {
        let mut symmetric = matrix;
        for (i, row) in symmetric.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = 0.5 * matrix[i][j] + 0.5 * matrix[j][i];
            }
        }
        Self::new(symmetric)
    }
    pub fn diagonal(variances: Vector6) -> Result<Self> {
        let mut matrix = [[0.0; 6]; 6];
        for (i, value) in variances.iter().enumerate() {
            matrix[i][i] = *value;
        }
        Self::new(matrix)
    }
    pub fn isotropic(position_std_m: f64, velocity_std_m_s: f64) -> Result<Self> {
        finite(position_std_m, "position uncertainty")?;
        finite(velocity_std_m_s, "velocity uncertainty")?;
        if position_std_m < 0.0 || velocity_std_m_s < 0.0 {
            return Err(invalid("uncertainties must be nonnegative"));
        }
        Self::diagonal([
            position_std_m.powi(2),
            position_std_m.powi(2),
            position_std_m.powi(2),
            velocity_std_m_s.powi(2),
            velocity_std_m_s.powi(2),
            velocity_std_m_s.powi(2),
        ])
    }
    pub fn matrix(self) -> Matrix6 {
        self.matrix
    }
    pub fn standard_deviations(self) -> Vector6 {
        std::array::from_fn(|i| self.matrix[i][i].sqrt())
    }
    pub fn propagated(self, jacobian: Matrix6) -> Result<Self> {
        for row in jacobian {
            for value in row {
                finite(value, "state Jacobian")?;
            }
        }
        Self::calculated(multiply(
            multiply(jacobian, self.matrix),
            transpose(jacobian),
        ))
    }
    /// Variance of a scalar linearized observable, h P hᵀ.
    pub fn projected_variance(self, jacobian: Vector6) -> Result<f64> {
        for value in jacobian {
            finite(value, "measurement Jacobian")?;
        }
        let ph = multiply_vector(self.matrix, jacobian);
        let variance: f64 = (0..6).map(|i| jacobian[i] * ph[i]).sum();
        finite(variance, "projected variance")?;
        if variance < 0.0 {
            return Err(invalid("projected covariance is negative"));
        }
        Ok(variance)
    }
    fn solve(self, rhs: Vector6) -> Result<Vector6> {
        let (normalized, scales) = correlation(self.matrix)?;
        if scales.contains(&0.0) {
            return Err(invalid(
                "smoothing requires nonsingular predicted covariance",
            ));
        }
        let l = cholesky(normalized, false)?;
        let mut y = [0.0; 6];
        for i in 0..6 {
            y[i] = (rhs[i] / scales[i] - (0..i).map(|j| l[i][j] * y[j]).sum::<f64>()) / l[i][i];
        }
        let mut x = [0.0; 6];
        for i in (0..6).rev() {
            x[i] = (y[i] - ((i + 1)..6).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
        }
        for i in 0..6 {
            x[i] /= scales[i];
            finite(x[i], "covariance solve")?;
        }
        Ok(x)
    }
}

/// Full state Jacobian of a deterministic rigid pose, including spin-induced
/// position/velocity covariance. Translation adds no uncertainty by itself.
pub fn pose_jacobian(pose: Pose) -> Result<Matrix6> {
    pose.angular_velocity_rad_s.validate()?;
    let basis = [
        Vec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
        Vec3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        },
        Vec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        },
    ];
    let mut j = [[0.0; 6]; 6];
    for (col, unit) in basis.into_iter().enumerate() {
        let r = pose.rotation.apply(unit);
        let spin = pose.angular_velocity_rad_s.cross(r);
        for (row, value) in [r.x, r.y, r.z].into_iter().enumerate() {
            j[row][col] = value;
            j[row + 3][col + 3] = value;
        }
        for (row, value) in [spin.x, spin.y, spin.z].into_iter().enumerate() {
            j[row + 3][col] = value;
        }
    }
    Ok(j)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrackState {
    pub epoch: Epoch,
    pub mean: State,
    pub covariance: Covariance6,
}
impl TrackState {
    pub fn new(epoch: Epoch, mean: State, covariance: Covariance6) -> Result<Self> {
        mean.position_m.validate()?;
        mean.velocity_m_s.validate()?;
        Ok(Self {
            epoch,
            mean,
            covariance,
        })
    }
    /// Constant-velocity prediction driven by continuous white acceleration.
    /// The three q values are acceleration-noise spectral densities, m²/s³,
    /// NOT acceleration standard deviations. q=0 gives a deterministic transition.
    pub fn predict(self, epoch: Epoch, acceleration_spectral_density: Vec3) -> Result<Self> {
        self.mean.position_m.validate()?;
        self.mean.velocity_m_s.validate()?;
        acceleration_spectral_density.validate()?;
        let q = [
            acceleration_spectral_density.x,
            acceleration_spectral_density.y,
            acceleration_spectral_density.z,
        ];
        if q.iter().any(|value| *value < 0.0) {
            return Err(invalid("process spectral densities must be nonnegative"));
        }
        let dt = epoch.duration_since(self.epoch);
        if dt < 0.0 {
            return Err(invalid("recursive prediction cannot run backward"));
        }
        let f = transition(dt);
        let mut process = [[0.0; 6]; 6];
        for (i, value) in q.into_iter().enumerate() {
            process[i][i] = value * dt.powi(3) / 3.0;
            process[i][i + 3] = value * dt.powi(2) / 2.0;
            process[i + 3][i] = process[i][i + 3];
            process[i + 3][i + 3] = value * dt;
        }
        let covariance = Covariance6::calculated(add(
            multiply(multiply(f, self.covariance.matrix), transpose(f)),
            process,
        ))?;
        Self::new(
            epoch,
            vector_state(multiply_vector(f, state_vector(self.mean)))?,
            covariance,
        )
    }
    pub fn transform(self, graph: &FrameGraph, from: FrameId, to: FrameId) -> Result<Self> {
        let pose = graph.relative_pose(from, to, self.epoch)?;
        Self::new(
            self.epoch,
            graph.transform(self.mean, from, to, self.epoch)?,
            self.covariance.propagated(pose_jacobian(pose)?)?,
        )
    }
    /// One scalar EKF update with a Joseph-form covariance update. Gate is a
    /// squared normalized innovation threshold (e.g. 9 for a scalar 3-sigma gate).
    /// Rejected measurements leave the predicted state and covariance unchanged.
    pub fn update(
        self,
        measurement: ScalarMeasurement,
        gate: Option<f64>,
    ) -> Result<(Self, UpdateDiagnostics)> {
        self.mean.position_m.validate()?;
        self.mean.velocity_m_s.validate()?;
        measurement.validate()?;
        if let Some(gate) = gate {
            finite(gate, "innovation gate")?;
            if gate <= 0.0 {
                return Err(invalid("innovation gate must be positive"));
            }
        }
        let residual = measurement.observed - measurement.predicted;
        finite(residual, "innovation")?;
        let r = measurement.noise_std.powi(2);
        finite(r, "measurement variance")?;
        if r == 0.0 {
            return Err(invalid("measurement variance underflows"));
        }
        let h = measurement.jacobian;
        let ph = multiply_vector(self.covariance.matrix, h);
        let innovation_variance = self.covariance.projected_variance(h)? + r;
        finite(innovation_variance, "innovation variance")?;
        let normalized_innovation_squared = (residual / innovation_variance.sqrt()).powi(2);
        finite(normalized_innovation_squared, "normalized innovation")?;
        let accepted = gate.is_none_or(|threshold| normalized_innovation_squared <= threshold);
        let diagnostics = UpdateDiagnostics {
            accepted,
            residual,
            innovation_variance,
            normalized_innovation_squared,
        };
        if !accepted {
            return Ok((self, diagnostics));
        }
        let gain: Vector6 = std::array::from_fn(|i| ph[i] / innovation_variance);
        let x = state_vector(self.mean);
        let mean = vector_state(std::array::from_fn(|i| x[i] + gain[i] * residual))?;
        let mut a = identity();
        let mut krk = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] -= gain[i] * h[j];
                krk[i][j] = gain[i] * r * gain[j];
            }
        }
        let covariance = Covariance6::calculated(add(
            multiply(multiply(a, self.covariance.matrix), transpose(a)),
            krk,
        ))?;
        Ok((Self::new(self.epoch, mean, covariance)?, diagnostics))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ScalarMeasurement {
    pub observed: f64,
    pub predicted: f64,
    pub jacobian: Vector6,
    pub noise_std: f64,
}
impl ScalarMeasurement {
    fn validate(self) -> Result<()> {
        finite(self.observed, "observation")?;
        finite(self.predicted, "prediction")?;
        finite(self.noise_std, "measurement noise")?;
        for value in self.jacobian {
            finite(value, "measurement Jacobian")?;
        }
        if self.noise_std <= 0.0 {
            return Err(invalid("measurement uncertainty must be positive"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy)]
pub struct UpdateDiagnostics {
    pub accepted: bool,
    pub residual: f64,
    pub innovation_variance: f64,
    pub normalized_innovation_squared: f64,
}

pub trait MeasurementModel {
    fn linearize(&self, predicted: &TrackState) -> Result<ScalarMeasurement>;
}
/// Exact scalar linear observable of the six-state vector; useful for Cartesian
/// position/velocity and already-decorrelated observation combinations.
pub struct LinearMeasurement {
    pub coefficients: Vector6,
    pub observed: f64,
    pub noise_std: f64,
}
impl MeasurementModel for LinearMeasurement {
    fn linearize(&self, predicted: &TrackState) -> Result<ScalarMeasurement> {
        let x = state_vector(predicted.mean);
        let measurement = ScalarMeasurement {
            observed: self.observed,
            predicted: (0..6).map(|i| self.coefficients[i] * x[i]).sum(),
            jacobian: self.coefficients,
            noise_std: self.noise_std,
        };
        measurement.validate()?;
        Ok(measurement)
    }
}

/// Absolute bistatic timing measurement, using retarded paths for both legs.
/// Endpoint worldlines, track state, and epoch must share one inertial reference.
/// The target is locally constant velocity during each light-time solve.
/// Central-difference steps must be chosen above local f64/solver resolution.
pub struct BistaticDelayMeasurement<'a> {
    pub transmitter: &'a dyn Worldline,
    pub receiver: &'a dyn Worldline,
    pub emission_epoch: Epoch,
    pub observed_delay_s: f64,
    pub noise_std_s: f64,
    pub known_delay_offset_s: f64,
    pub position_step_m: f64,
    pub velocity_step_m_s: f64,
    pub tolerance_s: f64,
}
impl MeasurementModel for BistaticDelayMeasurement<'_> {
    fn linearize(&self, predicted: &TrackState) -> Result<ScalarMeasurement> {
        if predicted.epoch != self.emission_epoch {
            return Err(invalid(
                "timing measurement needs a track at its emission epoch",
            ));
        }
        for value in [
            self.known_delay_offset_s,
            self.position_step_m,
            self.velocity_step_m_s,
            self.tolerance_s,
        ] {
            finite(value, "timing model parameter")?;
        }
        if self.position_step_m <= 0.0 || self.velocity_step_m_s <= 0.0 || self.tolerance_s <= 0.0 {
            return Err(invalid(
                "timing derivative steps/tolerance must be positive",
            ));
        }
        let evaluate = |x: Vector6| -> Result<f64> {
            let target = ConstantAcceleration {
                reference_epoch: predicted.epoch,
                state: vector_state(x)?,
                acceleration_m_s2: Vec3::ZERO,
            };
            let delay = bistatic_path(
                self.transmitter,
                &target,
                self.receiver,
                self.emission_epoch,
                self.tolerance_s,
            )?
            .total_vacuum_delay_s()
                + self.known_delay_offset_s;
            finite(delay, "bistatic prediction")?;
            Ok(delay)
        };
        let x = state_vector(predicted.mean);
        let mean = evaluate(x)?;
        let mut h = [0.0; 6];
        for i in 0..6 {
            let step = if i < 3 {
                self.position_step_m
            } else {
                self.velocity_step_m_s
            };
            let mut plus = x;
            let mut minus = x;
            plus[i] += step;
            minus[i] -= step;
            if plus[i] == x[i] || minus[i] == x[i] {
                return Err(invalid(
                    "timing derivative step is below coordinate precision",
                ));
            }
            h[i] = (evaluate(plus)? - evaluate(minus)?) / (2.0 * step);
        }
        let measurement = ScalarMeasurement {
            observed: self.observed_delay_s,
            predicted: mean,
            jacobian: h,
            noise_std: self.noise_std_s,
        };
        measurement.validate()?;
        Ok(measurement)
    }
}

/// Sequential filter with explicit observation IDs and chronological updates.
/// Failed operations are transactional; no state, epoch, or ID is consumed.
/// IDs of rejected outliers are consumed to prevent repeated reuse of evidence.
pub struct RecursiveTracker {
    current: TrackState,
    process_noise: Vec3,
    observations: HashSet<String>,
}
#[derive(Debug, Clone, Copy)]
pub struct Assimilation {
    pub predicted: TrackState,
    pub filtered: TrackState,
    pub diagnostics: UpdateDiagnostics,
}
impl RecursiveTracker {
    pub fn new(initial: TrackState, acceleration_spectral_density: Vec3) -> Result<Self> {
        initial.predict(initial.epoch, acceleration_spectral_density)?;
        Ok(Self {
            current: initial,
            process_noise: acceleration_spectral_density,
            observations: HashSet::new(),
        })
    }
    pub fn current(&self) -> TrackState {
        self.current
    }
    pub fn assimilate(
        &mut self,
        observation_id: impl Into<String>,
        epoch: Epoch,
        model: &dyn MeasurementModel,
        gate: Option<f64>,
    ) -> Result<Assimilation> {
        let id = observation_id.into();
        if id.trim().is_empty() || self.observations.contains(&id) {
            return Err(invalid("observation ID must be unique and nonempty"));
        }
        let predicted = self.current.predict(epoch, self.process_noise)?;
        let measurement = model.linearize(&predicted)?;
        let (filtered, diagnostics) = predicted.update(measurement, gate)?;
        self.current = filtered;
        self.observations.insert(id);
        Ok(Assimilation {
            predicted,
            filtered,
            diagnostics,
        })
    }
}

/// Rauch-Tung-Striebel fixed-interval smoother for the same constant-velocity
/// model used by [`TrackState::predict`]. `predictions[i]` is the prior at `filtered[i+1]`
/// obtained from `filtered[i]` BEFORE assimilating measurements at that epoch.
/// Store one final filtered state per distinct epoch, after all same-epoch updates.
/// This is Gaussian smoothing, not reprocessing or duplicating the observations.
pub fn smooth_history(
    filtered: &[TrackState],
    predictions: &[TrackState],
) -> Result<Vec<TrackState>> {
    if filtered.is_empty() {
        return Err(Error::NoData);
    }
    if predictions.len() != filtered.len() - 1 {
        return Err(invalid("smoother prediction/history lengths do not match"));
    }
    for state in filtered.iter().chain(predictions) {
        state.mean.position_m.validate()?;
        state.mean.velocity_m_s.validate()?;
    }
    for i in 0..predictions.len() {
        if filtered[i + 1].epoch <= filtered[i].epoch
            || predictions[i].epoch != filtered[i + 1].epoch
        {
            return Err(invalid(
                "smoother epochs must strictly increase and predictions must match",
            ));
        }
        // Validate a constant-velocity predicted mean; covariance may include process noise.
        let expected = multiply_vector(
            transition(filtered[i + 1].epoch.duration_since(filtered[i].epoch)),
            state_vector(filtered[i].mean),
        );
        let actual = state_vector(predictions[i].mean);
        for k in 0..6 {
            if (expected[k] - actual[k]).abs() > 1e-10 * expected[k].abs().max(1.0) {
                return Err(invalid(
                    "smoother prediction must use the constant-velocity transition",
                ));
            }
        }
    }
    let mut smoothed = filtered.to_vec();
    for i in (0..predictions.len()).rev() {
        let dt = filtered[i + 1].epoch.duration_since(filtered[i].epoch);
        let cross = multiply(filtered[i].covariance.matrix, transpose(transition(dt)));
        let mut gain = [[0.0; 6]; 6];
        for row in 0..6 {
            gain[row] = predictions[i].covariance.solve(cross[row])?;
        }
        let next = state_vector(smoothed[i + 1].mean);
        let predicted = state_vector(predictions[i].mean);
        let residual = std::array::from_fn(|k| next[k] - predicted[k]);
        let correction = multiply_vector(gain, residual);
        let current = state_vector(filtered[i].mean);
        let covariance = Covariance6::calculated(add(
            filtered[i].covariance.matrix,
            multiply(
                multiply(
                    gain,
                    subtract(
                        smoothed[i + 1].covariance.matrix,
                        predictions[i].covariance.matrix,
                    ),
                ),
                transpose(gain),
            ),
        ))?;
        smoothed[i] = TrackState::new(
            filtered[i].epoch,
            vector_state(std::array::from_fn(|k| current[k] + correction[k]))?,
            covariance,
        )?;
    }
    Ok(smoothed)
}

/// C1 interpolation of smoothed means. Covariances remain at their knots; this
/// does not invent an uncertainty interpolation or change the process model.
pub fn mean_worldpath(history: &[TrackState]) -> Result<HermiteWorldline> {
    HermiteWorldline::new(
        history
            .iter()
            .map(|state| (state.epoch, state.mean))
            .collect(),
    )
}
