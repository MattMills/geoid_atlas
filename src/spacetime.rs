//! Continuous worldlines, retarded bistatic paths, and phase-aligned integration.
//!
//! These are forward-model primitives, not an automatic orbit-determination or
//! sensor-fusion engine. All paths share one declared coordinate time and frame.
use crate::frames::{Epoch, FrameGraph, FrameId, State, Vec3};
use crate::interferometry::SPEED_OF_LIGHT_M_S;
use crate::{Error, Result, finite};
use std::collections::{BTreeMap, HashSet};

pub trait Worldline: Send + Sync {
    fn state_at(&self, epoch: Epoch) -> Result<State>;
}

#[derive(Debug, Clone, Copy)]
pub struct ConstantAcceleration {
    pub reference_epoch: Epoch,
    pub state: State,
    pub acceleration_m_s2: Vec3,
}
impl Worldline for ConstantAcceleration {
    fn state_at(&self, epoch: Epoch) -> Result<State> {
        self.state.position_m.validate()?;
        self.state.velocity_m_s.validate()?;
        self.acceleration_m_s2.validate()?;
        let dt = epoch.duration_since(self.reference_epoch);
        let state = State {
            position_m: self.state.position_m
                + self.state.velocity_m_s * dt
                + self.acceleration_m_s2 * (0.5 * dt * dt),
            velocity_m_s: self.state.velocity_m_s + self.acceleration_m_s2 * dt,
        };
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        Ok(state)
    }
}

/// Apply the graph at every queried epoch, including each light-time iteration.
pub struct FramedWorldline<'a> {
    pub graph: &'a FrameGraph,
    pub frame: FrameId,
    pub local: &'a dyn Worldline,
}
impl Worldline for FramedWorldline<'_> {
    fn state_at(&self, epoch: Epoch) -> Result<State> {
        self.graph.transform(
            self.local.state_at(epoch)?,
            self.frame,
            self.graph.root(),
            epoch,
        )
    }
}

/// Cubic Hermite worldpath through timestamped position/velocity knots.
/// Position and first derivative are contiguous at each knot; acceleration may jump.
/// No extrapolation: use a process model or append new knots instead.
pub struct HermiteWorldline {
    knots: Vec<(Epoch, State)>,
}
impl HermiteWorldline {
    pub fn new(knots: Vec<(Epoch, State)>) -> Result<Self> {
        if knots.len() < 2 {
            return Err(Error::InvalidInput(
                "worldpath needs at least two knots".into(),
            ));
        }
        for (i, (epoch, state)) in knots.iter().enumerate() {
            state.position_m.validate()?;
            state.velocity_m_s.validate()?;
            if i > 0 && *epoch <= knots[i - 1].0 {
                return Err(Error::InvalidInput(
                    "worldpath epochs must strictly increase".into(),
                ));
            }
        }
        Ok(Self { knots })
    }
}
impl Worldline for HermiteWorldline {
    fn state_at(&self, epoch: Epoch) -> Result<State> {
        if epoch < self.knots[0].0 || epoch > self.knots[self.knots.len() - 1].0 {
            return Err(Error::OutsideCoverage);
        }
        let i = self
            .knots
            .partition_point(|(t, _)| *t <= epoch)
            .saturating_sub(1)
            .min(self.knots.len() - 2);
        let (t0, s0) = self.knots[i];
        let (t1, s1) = self.knots[i + 1];
        let dt = t1.duration_since(t0);
        let u = epoch.duration_since(t0) / dt;
        let position = s0.position_m * (2.0 * u.powi(3) - 3.0 * u * u + 1.0)
            + s0.velocity_m_s * (dt * (u.powi(3) - 2.0 * u * u + u))
            + s1.position_m * (-2.0 * u.powi(3) + 3.0 * u * u)
            + s1.velocity_m_s * (dt * (u.powi(3) - u * u));
        let velocity = s0.position_m * ((6.0 * u * u - 6.0 * u) / dt)
            + s0.velocity_m_s * (3.0 * u * u - 4.0 * u + 1.0)
            + s1.position_m * ((-6.0 * u * u + 6.0 * u) / dt)
            + s1.velocity_m_s * (3.0 * u * u - 2.0 * u);
        position.validate()?;
        velocity.validate()?;
        Ok(State {
            position_m: position,
            velocity_m_s: velocity,
        })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RetardedEvent {
    pub epoch: Epoch,
    pub state: State,
    pub travel_time_s: f64,
    pub iterations: u32,
}

/// Vacuum light travel to a moving worldline, for one fixed emission event.
pub fn retarded_event(
    emitter_position_m: Vec3,
    emission_epoch: Epoch,
    receiver: &dyn Worldline,
    tolerance_s: f64,
) -> Result<RetardedEvent> {
    emitter_position_m.validate()?;
    finite(tolerance_s, "light-time tolerance")?;
    if tolerance_s <= 0.0 {
        return Err(Error::InvalidInput("tolerance must be positive".into()));
    }
    let mut dt = 0.0;
    for iterations in 1..=200 {
        let event_epoch = emission_epoch.shifted(dt)?;
        let state = receiver.state_at(event_epoch)?;
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        if state.velocity_m_s.norm() >= SPEED_OF_LIGHT_M_S {
            return Err(Error::InvalidInput("worldline must be subluminal".into()));
        }
        let next = (state.position_m - emitter_position_m).norm() / SPEED_OF_LIGHT_M_S;
        finite(next, "light time")?;
        if (next - dt).abs() <= tolerance_s {
            let epoch = emission_epoch.shifted(next)?;
            let state = receiver.state_at(epoch)?;
            state.position_m.validate()?;
            state.velocity_m_s.validate()?;
            return Ok(RetardedEvent {
                epoch,
                state,
                travel_time_s: next,
                iterations,
            });
        }
        dt = next;
    }
    Err(Error::InvalidInput(
        "retarded worldline solve did not converge".into(),
    ))
}

#[derive(Debug, Clone, Copy)]
pub struct BistaticPath {
    pub emission_epoch: Epoch,
    pub transmitter: State,
    pub scattering: RetardedEvent,
    pub reception: RetardedEvent,
}
impl BistaticPath {
    pub fn total_vacuum_delay_s(self) -> f64 {
        self.scattering.travel_time_s + self.reception.travel_time_s
    }
}
/// Transmitter -> moving target -> moving receiver, solving both legs at their
/// distinct epochs. Reflectivity, medium, occultation, and clocks are separate models.
pub fn bistatic_path(
    transmitter: &dyn Worldline,
    target: &dyn Worldline,
    receiver: &dyn Worldline,
    emission_epoch: Epoch,
    tolerance_s: f64,
) -> Result<BistaticPath> {
    let tx = transmitter.state_at(emission_epoch)?;
    tx.position_m.validate()?;
    tx.velocity_m_s.validate()?;
    if tx.velocity_m_s.norm() >= SPEED_OF_LIGHT_M_S {
        return Err(Error::InvalidInput("transmitter must be subluminal".into()));
    }
    let scattering = retarded_event(tx.position_m, emission_epoch, target, tolerance_s)?;
    let reception = retarded_event(
        scattering.state.position_m,
        scattering.epoch,
        receiver,
        tolerance_s,
    )?;
    Ok(BistaticPath {
        emission_epoch,
        transmitter: tx,
        scattering,
        reception,
    })
}

/// Calibrated complex baseband sample within ONE phase-linked channel.
/// phase_rad is the forward-model phase relative to the acquisition reference,
/// not an absolute astronomical carrier phase computed from a large Julian date.
#[derive(Debug, Clone, Copy)]
pub struct PhaseSample {
    pub real: f64,
    pub imaginary: f64,
    pub predicted_phase_rad: f64,
    /// Independent Gaussian noise standard deviation in each quadrature.
    pub noise_std: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct CoherentEstimate {
    pub real: f64,
    pub imaginary: f64,
    pub noise_std: f64,
    pub samples: usize,
}

/// Phase-align and inverse-variance combine independent samples of a constant
/// complex reflectivity. Use separate calls per phase-linked channel/modality.
/// Time-varying reflectivity requires a shorter local window or a different model.
pub fn integrate_coherent(samples: &[PhaseSample]) -> Result<CoherentEstimate> {
    if samples.is_empty() {
        return Err(Error::NoData);
    }
    let mut min_sigma = f64::INFINITY;
    for s in samples {
        for v in [s.real, s.imaginary, s.predicted_phase_rad, s.noise_std] {
            finite(v, "coherent sample")?;
        }
        if s.noise_std <= 0.0 {
            return Err(Error::InvalidInput("sample noise must be positive".into()));
        }
        min_sigma = min_sigma.min(s.noise_std);
    }
    let weight_sum: f64 = samples
        .iter()
        .map(|s| (min_sigma / s.noise_std).powi(2))
        .sum();
    let (mut real, mut imaginary) = (0.0, 0.0);
    for s in samples {
        let (sin, cos) = s.predicted_phase_rad.sin_cos();
        let weight = (min_sigma / s.noise_std).powi(2) / weight_sum;
        real += (s.real * cos + s.imaginary * sin) * weight;
        imaginary += (s.imaginary * cos - s.real * sin) * weight;
    }
    finite(real, "integrated real")?;
    finite(imaginary, "integrated imaginary")?;
    Ok(CoherentEstimate {
        real,
        imaginary,
        noise_std: min_sigma / weight_sum.sqrt(),
        samples: samples.len(),
    })
}

/// Local hypothesis search cell. Refinement is explicit and bounded: halve
/// spatial/temporal extents about the same centre, without constructing a global volume.
#[derive(Debug, Clone, Copy)]
pub struct LocalCell {
    centre_m: Vec3,
    centre_epoch: Epoch,
    half_size_m: f64,
    half_duration_s: f64,
}
impl LocalCell {
    pub fn new(
        centre_m: Vec3,
        centre_epoch: Epoch,
        half_size_m: f64,
        half_duration_s: f64,
    ) -> Result<Self> {
        centre_m.validate()?;
        finite(half_size_m, "cell half size")?;
        finite(half_duration_s, "cell duration")?;
        if half_size_m <= 0.0 || half_duration_s <= 0.0 {
            return Err(Error::InvalidInput("cell extents must be positive".into()));
        }
        Ok(Self {
            centre_m,
            centre_epoch,
            half_size_m,
            half_duration_s,
        })
    }
    pub fn contains(self, position_m: Vec3, epoch: Epoch) -> Result<bool> {
        position_m.validate()?;
        let d = position_m - self.centre_m;
        Ok(d.x.abs() <= self.half_size_m
            && d.y.abs() <= self.half_size_m
            && d.z.abs() <= self.half_size_m
            && epoch.duration_since(self.centre_epoch).abs() <= self.half_duration_s)
    }
    pub fn refined(self) -> Result<Self> {
        Self::new(
            self.centre_m,
            self.centre_epoch,
            self.half_size_m / 2.0,
            self.half_duration_s / 2.0,
        )
    }
    pub fn half_size_m(self) -> f64 {
        self.half_size_m
    }
    pub fn half_duration_s(self) -> f64 {
        self.half_duration_s
    }
}

/// Sensor-specific likelihood in a declared common coordinate system.
/// Optical astrometry, photon counts, weather and complex RF require DIFFERENT
/// likelihood implementations, even when they constrain the same worldpath.
pub trait NegativeLogLikelihood: Send + Sync {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64>;
}

/// Independent-axis Gaussian state prior at one epoch in the same frame as the
/// worldpath. Can carry a previous window's estimate with process-noise inflation.
/// Cross-axis and position/velocity covariance require a custom joint model.
pub struct GaussianStatePrior {
    pub epoch: Epoch,
    pub mean: State,
    pub position_std_m: Vec3,
    pub velocity_std_m_s: Vec3,
}
impl NegativeLogLikelihood for GaussianStatePrior {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64> {
        let state = trajectory.state_at(self.epoch)?;
        self.mean.position_m.validate()?;
        self.mean.velocity_m_s.validate()?;
        let dp = state.position_m - self.mean.position_m;
        let dv = state.velocity_m_s - self.mean.velocity_m_s;
        let mut cost = 0.0;
        for (delta, sigma) in [
            (dp.x, self.position_std_m.x),
            (dp.y, self.position_std_m.y),
            (dp.z, self.position_std_m.z),
            (dv.x, self.velocity_std_m_s.x),
            (dv.y, self.velocity_std_m_s.y),
            (dv.z, self.velocity_std_m_s.z),
        ] {
            finite(delta, "state residual")?;
            finite(sigma, "prior uncertainty")?;
            if sigma <= 0.0 {
                return Err(Error::InvalidInput(
                    "prior uncertainty must be positive".into(),
                ));
            }
            cost += 0.5 * (delta / sigma).powi(2);
        }
        finite(cost, "state prior cost")?;
        Ok(cost)
    }
}

/// Poisson negative log likelihood without the observation-only log(n!) term.
/// Used for photon-count observations (e.g. X-rays), not complex field coherence.
pub fn poisson_count_cost(observed_count: u64, expected_count: f64) -> Result<f64> {
    finite(expected_count, "expected count")?;
    if expected_count < 0.0 {
        return Err(Error::InvalidInput(
            "expected count must be nonnegative".into(),
        ));
    }
    if expected_count == 0.0 {
        return if observed_count == 0 {
            Ok(0.0)
        } else {
            Err(Error::InvalidInput(
                "nonzero count has zero model probability".into(),
            ))
        };
    }
    let cost = expected_count - observed_count as f64 * expected_count.ln();
    finite(cost, "count likelihood")?;
    Ok(cost)
}

struct LikelihoodTerm<'a> {
    id: String,
    group: String,
    model: Box<dyn NegativeLogLikelihood + 'a>,
}
#[derive(Default)]
pub struct JointLikelihood<'a> {
    terms: Vec<LikelihoodTerm<'a>>,
}
impl<'a> JointLikelihood<'a> {
    /// Group names declare conditional independence of likelihood blocks.
    /// Shared calibrations must be modeled jointly within a single block.
    pub fn add(
        &mut self,
        id: impl Into<String>,
        independence_group: impl Into<String>,
        model: impl NegativeLogLikelihood + 'a,
    ) -> Result<()> {
        let id = id.into();
        let group = independence_group.into();
        if id.trim().is_empty()
            || group.trim().is_empty()
            || self.terms.iter().any(|t| t.id == id || t.group == group)
        {
            return Err(Error::InvalidInput(
                "likelihood IDs and independence groups must be unique and nonempty".into(),
            ));
        }
        self.terms.push(LikelihoodTerm {
            id,
            group,
            model: Box::new(model),
        });
        Ok(())
    }
}
impl NegativeLogLikelihood for JointLikelihood<'_> {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64> {
        if self.terms.is_empty() {
            return Err(Error::NoData);
        }
        let mut total = 0.0;
        for term in &self.terms {
            let cost = term.model.evaluate(trajectory)?;
            finite(cost, "likelihood cost")?;
            total += cost;
        }
        finite(total, "joint likelihood")?;
        Ok(total)
    }
}

#[derive(Debug, Clone)]
pub struct BistaticObservation {
    pub id: String,
    pub coherent_channel: String,
    pub emission_epoch: Epoch,
    pub transmitter_index: usize,
    pub receiver_index: usize,
    pub frequency_hz: f64,
    /// Reference-model vacuum travel time used by the acquisition/downconversion.
    pub reference_delay_s: f64,
    pub real: f64,
    pub imaginary: f64,
    pub noise_std: f64,
}

/// Profile likelihood of constant complex reflectivity per calibrated coherent
/// channel, using independent quadrature noise. Separate channel names prevent
/// RF/optical or uncalibrated instruments from being summed as coherent fields.
pub struct CoherentBistaticLikelihood<'a> {
    transmitters: Vec<&'a dyn Worldline>,
    receivers: Vec<&'a dyn Worldline>,
    observations: Vec<BistaticObservation>,
    tolerance_s: f64,
}
impl<'a> CoherentBistaticLikelihood<'a> {
    pub fn new(
        transmitters: Vec<&'a dyn Worldline>,
        receivers: Vec<&'a dyn Worldline>,
        observations: Vec<BistaticObservation>,
        tolerance_s: f64,
    ) -> Result<Self> {
        finite(tolerance_s, "tolerance")?;
        if tolerance_s <= 0.0 || observations.is_empty() {
            return Err(Error::InvalidInput(
                "require observations and positive tolerance".into(),
            ));
        }
        let mut ids = HashSet::new();
        for o in &observations {
            if o.id.trim().is_empty()
                || o.coherent_channel.trim().is_empty()
                || !ids.insert(o.id.clone())
                || o.transmitter_index >= transmitters.len()
                || o.receiver_index >= receivers.len()
            {
                return Err(Error::InvalidInput(
                    "observation IDs/channels/endpoint indices are invalid".into(),
                ));
            }
            for v in [
                o.frequency_hz,
                o.reference_delay_s,
                o.real,
                o.imaginary,
                o.noise_std,
            ] {
                finite(v, "observation parameter")?;
            }
            if o.frequency_hz <= 0.0 || o.reference_delay_s < 0.0 || o.noise_std <= 0.0 {
                return Err(Error::InvalidInput(
                    "require f>0, reference delay>=0, noise>0".into(),
                ));
            }
        }
        Ok(Self {
            transmitters,
            receivers,
            observations,
            tolerance_s,
        })
    }
}
impl NegativeLogLikelihood for CoherentBistaticLikelihood<'_> {
    fn evaluate(&self, trajectory: &dyn Worldline) -> Result<f64> {
        let mut channels: BTreeMap<&str, Vec<PhaseSample>> = BTreeMap::new();
        for observation in &self.observations {
            let path = bistatic_path(
                self.transmitters[observation.transmitter_index],
                trajectory,
                self.receivers[observation.receiver_index],
                observation.emission_epoch,
                self.tolerance_s,
            )?;
            let residual_delay = path.total_vacuum_delay_s() - observation.reference_delay_s;
            let phase = 2.0
                * std::f64::consts::PI
                * (observation.frequency_hz * residual_delay).rem_euclid(1.0);
            finite(phase, "residual phase")?;
            channels
                .entry(&observation.coherent_channel)
                .or_default()
                .push(PhaseSample {
                    real: observation.real,
                    imaginary: observation.imaginary,
                    predicted_phase_rad: phase,
                    noise_std: observation.noise_std,
                });
        }
        let mut cost = 0.0;
        for samples in channels.values() {
            let amplitude = integrate_coherent(samples)?;
            for s in samples {
                let (sin, cos) = s.predicted_phase_rad.sin_cos();
                let dr = (s.real * cos + s.imaginary * sin - amplitude.real) / s.noise_std;
                let di = (s.imaginary * cos - s.real * sin - amplitude.imaginary) / s.noise_std;
                cost += 0.5 * (dr * dr + di * di);
            }
        }
        finite(cost, "coherent likelihood")?;
        Ok(cost)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LocalSearch {
    pub position_half_width_m: f64,
    pub velocity_half_width_m_s: f64,
    pub levels: u32,
    pub sweeps_per_level: u32,
}
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFit {
    pub trajectory: ConstantAcceleration,
    pub negative_log_likelihood: f64,
    pub evaluations: usize,
}

/// Bounded coarse-to-fine coordinate search over position and velocity (six
/// parameters). Acceleration is fixed. This is a local fit, not a global optimum,
/// posterior covariance, recursive Bayesian filter, or phase-ambiguity resolver.
/// Warm-start with the previous window's fit for successive local windows.
pub fn refine_trajectory(
    initial: ConstantAcceleration,
    settings: LocalSearch,
    likelihood: &dyn NegativeLogLikelihood,
) -> Result<TrajectoryFit> {
    for v in [
        settings.position_half_width_m,
        settings.velocity_half_width_m_s,
    ] {
        finite(v, "search width")?;
    }
    if settings.position_half_width_m <= 0.0
        || settings.velocity_half_width_m_s <= 0.0
        || settings.levels == 0
        || settings.levels > 40
        || settings.sweeps_per_level == 0
        || settings.sweeps_per_level > 100
    {
        return Err(Error::InvalidInput(
            "invalid bounded search settings".into(),
        ));
    }
    initial.state_at(initial.reference_epoch)?;
    let mut parameters = [
        initial.state.position_m.x,
        initial.state.position_m.y,
        initial.state.position_m.z,
        initial.state.velocity_m_s.x,
        initial.state.velocity_m_s.y,
        initial.state.velocity_m_s.z,
    ];
    let origin = parameters;
    let widths = [
        settings.position_half_width_m,
        settings.position_half_width_m,
        settings.position_half_width_m,
        settings.velocity_half_width_m_s,
        settings.velocity_half_width_m_s,
        settings.velocity_half_width_m_s,
    ];
    let make = |p: [f64; 6]| ConstantAcceleration {
        state: State {
            position_m: Vec3 {
                x: p[0],
                y: p[1],
                z: p[2],
            },
            velocity_m_s: Vec3 {
                x: p[3],
                y: p[4],
                z: p[5],
            },
        },
        ..initial
    };
    let mut best = likelihood.evaluate(&initial)?;
    finite(best, "initial likelihood")?;
    let mut evaluations = 1;
    for level in 0..settings.levels {
        for _ in 0..settings.sweeps_per_level {
            let mut changed = false;
            for axis in 0..6 {
                let starting = parameters;
                let step = widths[axis] * 2.0_f64.powi(-(level as i32));
                for sign in [-1.0, 1.0] {
                    let mut candidate = starting;
                    candidate[axis] += sign * step;
                    if (candidate[axis] - origin[axis]).abs() > widths[axis] * (1.0 + 1e-12) {
                        continue;
                    }
                    let value = likelihood.evaluate(&make(candidate))?;
                    finite(value, "candidate likelihood")?;
                    evaluations += 1;
                    if value < best {
                        best = value;
                        parameters = candidate;
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }
    Ok(TrajectoryFit {
        trajectory: make(parameters),
        negative_log_likelihood: best,
        evaluations,
    })
}
