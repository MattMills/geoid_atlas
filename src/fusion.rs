//! Sensor-specific forward predictions and correlated multimodal likelihoods.
//! Correlation is statistical agreement with a hypothesis, not cross-band field
//! phase coherence. Environmental nuisance parameters may be captured by providers;
//! to jointly estimate them, supply a larger hypothesis/optimizer outside Worldline.
use crate::frames::{Epoch, Vec3};
use crate::spacetime::{NegativeLogLikelihood, Worldline, poisson_count_cost};
use crate::{Error, Result, finite};

/// Predict a block of observables in the EXACT units/order of its observations.
/// Implementations carry their timestamps, frame adapters, calibration, coverage,
/// and response models. A closure can capture an atlas or an environmental model.
pub trait ObservableModel: Send + Sync {
    fn predict(&self, trajectory: &dyn Worldline) -> Result<Vec<f64>>;
}
impl<F> ObservableModel for F
where
    F: Fn(&dyn Worldline) -> Result<Vec<f64>> + Send + Sync,
{
    fn predict(&self, trajectory: &dyn Worldline) -> Result<Vec<f64>> {
        self(trajectory)
    }
}

/// Full correlated Gaussian block, including cross-sensor covariance. Uses
/// correlation-space Cholesky so e.g. metres and seconds can coexist numerically.
/// Cost omits constants for FIXED covariance; hypothesis-dependent covariance
/// requires its log determinant in a custom likelihood. Singular blocks fail.
pub struct GaussianObservation<M> {
    observed: Vec<f64>,
    scales: Vec<f64>,
    cholesky: Vec<Vec<f64>>,
    model: M,
}
impl<M: ObservableModel> GaussianObservation<M> {
    pub fn new(observed: Vec<f64>, covariance: Vec<Vec<f64>>, model: M) -> Result<Self> {
        let n = observed.len();
        if n == 0 || covariance.len() != n || covariance.iter().any(|row| row.len() != n) {
            return Err(Error::InvalidInput(
                "nonempty matching observation/covariance dimensions required".into(),
            ));
        }
        for &value in &observed {
            finite(value, "observed value")?;
        }
        let mut scales = Vec::with_capacity(n);
        for (i, row) in covariance.iter().enumerate() {
            finite(row[i], "observation variance")?;
            if row[i] <= 0.0 {
                return Err(Error::InvalidInput(
                    "observation variances must be positive".into(),
                ));
            }
            scales.push(row[i].sqrt());
        }
        let mut correlation = vec![vec![0.0; n]; n];
        for (i, row) in correlation.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                finite(covariance[i][j], "observation covariance")?;
                let a = covariance[i][j] / scales[i] / scales[j];
                let b = covariance[j][i] / scales[j] / scales[i];
                finite(a, "observation correlation")?;
                finite(b, "observation correlation")?;
                if (a - b).abs() > 1e-12 {
                    return Err(Error::InvalidInput(
                        "observation covariance must be symmetric".into(),
                    ));
                }
                *value = 0.5 * (a + b);
            }
        }
        let mut cholesky = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..=i {
                let residual = correlation[i][j]
                    - (0..j).map(|k| cholesky[i][k] * cholesky[j][k]).sum::<f64>();
                finite(residual, "observation factorization")?;
                if i == j {
                    if residual <= 1e-14 {
                        return Err(Error::InvalidInput(
                            "observation covariance must be positive definite and well conditioned"
                                .into(),
                        ));
                    }
                    cholesky[i][j] = residual.sqrt();
                } else {
                    cholesky[i][j] = residual / cholesky[j][j];
                }
            }
        }
        Ok(Self {
            observed,
            scales,
            cholesky,
            model,
        })
    }
}
impl<M: ObservableModel> NegativeLogLikelihood for GaussianObservation<M> {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64> {
        let predicted = self.model.predict(trajectory)?;
        if predicted.len() != self.observed.len() {
            return Err(Error::InvalidInput("prediction dimensions changed".into()));
        }
        let mut whitened = Vec::with_capacity(predicted.len());
        for (i, &prediction) in predicted.iter().enumerate() {
            finite(prediction, "predicted observable")?;
            let residual = (self.observed[i] - prediction) / self.scales[i];
            let value = (residual
                - (0..i)
                    .map(|j| self.cholesky[i][j] * whitened[j])
                    .sum::<f64>())
                / self.cholesky[i][i];
            finite(value, "whitened observation residual")?;
            whitened.push(value);
        }
        let cost = 0.5 * whitened.iter().map(|v| v * v).sum::<f64>();
        finite(cost, "Gaussian observation cost")?;
        Ok(cost)
    }
}

/// Independent photon counts conditional on the response/exposure/background
/// provider. Correlated counts require a joint model, not this Poisson adapter.
pub struct PhotonObservation<M> {
    pub counts: Vec<u64>,
    pub model: M,
}
impl<M: ObservableModel> NegativeLogLikelihood for PhotonObservation<M> {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64> {
        let expected = self.model.predict(trajectory)?;
        if self.counts.is_empty() || expected.len() != self.counts.len() {
            return Err(Error::InvalidInput(
                "nonempty matching photon bins required".into(),
            ));
        }
        let mut cost = 0.0;
        for (&count, expectation) in self.counts.iter().zip(expected) {
            cost += poisson_count_cost(count, expectation)?;
        }
        finite(cost, "photon observation cost")?;
        Ok(cost)
    }
}

/// Backward light-time solve: target EMISSION event corresponding to a fixed
/// sensor reception event. Vacuum only; clocks/medium/aberration must be supplied.
pub fn emission_epoch(
    target: &dyn Worldline,
    receiver_m: Vec3,
    reception_epoch: Epoch,
    tolerance_s: f64,
) -> Result<Epoch> {
    receiver_m.validate()?;
    finite(tolerance_s, "light-time tolerance")?;
    if tolerance_s <= 0.0 {
        return Err(Error::InvalidInput(
            "light-time tolerance must be positive".into(),
        ));
    }
    let mut delay = 0.0;
    for _ in 0..200 {
        let epoch = reception_epoch.shifted(-delay)?;
        let state = target.state_at(epoch)?;
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        if state.velocity_m_s.norm() >= crate::interferometry::SPEED_OF_LIGHT_M_S {
            return Err(Error::InvalidInput("worldline must be subluminal".into()));
        }
        let next =
            (state.position_m - receiver_m).norm() / crate::interferometry::SPEED_OF_LIGHT_M_S;
        finite(next, "light time")?;
        if (next - delay).abs() <= tolerance_s {
            return reception_epoch.shifted(-next);
        }
        delay = next;
    }
    Err(Error::InvalidInput(
        "backward light-time solve did not converge".into(),
    ))
}

/// Image/astrometric bearing likelihood in the observed direction's tangent plane.
/// Small-angle, independent isotropic uncertainty. Intrinsics/distortion/attitude
/// must first convert a pixel to a calibrated common-frame line of sight. This
/// model samples target emission time, sensor reception time; no relativistic aberration.
pub struct BearingObservation<'a> {
    pub sensor: &'a dyn Worldline,
    pub reception_epoch: Epoch,
    pub observed_direction: Vec3,
    pub angular_std_rad: f64,
    pub light_time_tolerance_s: f64,
}
impl NegativeLogLikelihood for BearingObservation<'_> {
    fn evaluate(&self, target: &dyn Worldline) -> Result<f64> {
        finite(self.angular_std_rad, "bearing uncertainty")?;
        if self.angular_std_rad <= 0.0 {
            return Err(Error::InvalidInput(
                "bearing uncertainty must be positive".into(),
            ));
        }
        let observed = self.observed_direction.normalized()?;
        let sensor = self.sensor.state_at(self.reception_epoch)?;
        let emitted = emission_epoch(
            target,
            sensor.position_m,
            self.reception_epoch,
            self.light_time_tolerance_s,
        )?;
        let predicted = (target.state_at(emitted)?.position_m - sensor.position_m).normalized()?;
        let forward = predicted.dot(observed);
        if forward <= 0.0 {
            return Err(Error::InvalidInput(
                "predicted bearing lies behind tangent plane".into(),
            ));
        }
        let tangent = predicted - observed * forward;
        let cost = 0.5 * (tangent.norm() / forward / self.angular_std_rad).powi(2);
        finite(cost, "bearing cost")?;
        Ok(cost)
    }
}
