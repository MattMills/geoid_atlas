//! Signed delay convention: arrival at receiver 2 minus arrival at receiver 1.
//! Geometry is classical; propagation, clock, and gravitational corrections are explicit.
use crate::frames::{Epoch, FrameGraph, FrameId, State, Vec3};
use crate::{Error, Result, finite};

pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

#[derive(Debug, Clone, Copy)]
pub struct FarFieldDelay {
    pub seconds: f64,
    pub rate_s_per_s: f64,
    /// Baseline coordinates in metres; w points toward the source.
    pub uvw_m: Vec3,
}

/// Both states and the source direction must be expressed in the same inertial frame.
/// Source direction points FROM the array TOWARD the source.
/// Delay rate holds that direction fixed; source-direction motion requires an
/// additional -baseline dot direction_derivative / c term.
pub fn far_field(
    receiver1: State,
    receiver2: State,
    source_direction: Vec3,
) -> Result<FarFieldDelay> {
    receiver1.position_m.validate()?;
    receiver2.position_m.validate()?;
    receiver1.velocity_m_s.validate()?;
    receiver2.velocity_m_s.validate()?;
    let w = source_direction.normalized()?;
    let reference = if w.z.abs() > 0.99 {
        Vec3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    } else {
        Vec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    };
    let u = reference.cross(w).normalized()?;
    let v = w.cross(u);
    let baseline = receiver2.position_m - receiver1.position_m;
    Ok(FarFieldDelay {
        seconds: -baseline.dot(w) / SPEED_OF_LIGHT_M_S,
        rate_s_per_s: -(receiver2.velocity_m_s - receiver1.velocity_m_s).dot(w)
            / SPEED_OF_LIGHT_M_S,
        uvw_m: Vec3 {
            x: baseline.dot(u),
            y: baseline.dot(v),
            z: baseline.dot(w),
        },
    })
}

/// Finite-distance spherical wavefront at fixed receiver positions.
/// Uses a difference-of-squares expression to reduce cancellation for distant sources.
pub fn near_field_static(source: Vec3, receiver1: Vec3, receiver2: Vec3) -> Result<f64> {
    source.validate()?;
    receiver1.validate()?;
    receiver2.validate()?;
    let a = receiver1 - source;
    let b = receiver2 - source;
    let total = a.norm() + b.norm();
    if total == 0.0 {
        return Ok(0.0);
    }
    let baseline = receiver2 - receiver1;
    let delay = baseline.dot(a + b) / total / SPEED_OF_LIGHT_M_S;
    finite(delay, "delay")?;
    Ok(delay)
}

/// Receiver state at reference_epoch in its own frame. Body rotation/orbit comes
/// from the graph; local motion is linearly extrapolated during light travel.
#[derive(Debug, Clone, Copy)]
pub struct Receiver {
    pub frame: FrameId,
    pub reference_epoch: Epoch,
    pub local_state: State,
}
impl Receiver {
    pub fn state_in_root(self, graph: &FrameGraph, epoch: Epoch) -> Result<State> {
        let dt = epoch.duration_since(self.reference_epoch);
        graph.transform(
            State {
                position_m: self.local_state.position_m + self.local_state.velocity_m_s * dt,
                velocity_m_s: self.local_state.velocity_m_s,
            },
            self.frame,
            graph.root(),
            epoch,
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Arrival {
    pub epoch: Epoch,
    pub travel_time_s: f64,
    pub iterations: u32,
}

/// Solve vacuum light time for a specified emission event and a moving receiver.
/// The source position is in the common root at emission. No extrapolation of
/// the source after emission is needed for that wavefront. Requires subluminal motion.
pub fn arrival(
    graph: &FrameGraph,
    source_at_emission: Vec3,
    emission: Epoch,
    receiver: Receiver,
    tolerance_s: f64,
) -> Result<Arrival> {
    source_at_emission.validate()?;
    finite(tolerance_s, "tolerance")?;
    if tolerance_s <= 0.0 {
        return Err(Error::InvalidInput("tolerance must be positive".into()));
    }
    let mut dt = 0.0;
    for iterations in 1..=200 {
        let state = receiver.state_in_root(graph, emission.shifted(dt)?)?;
        if state.velocity_m_s.norm() >= SPEED_OF_LIGHT_M_S {
            return Err(Error::InvalidInput("receiver must be subluminal".into()));
        }
        let next = (state.position_m - source_at_emission).norm() / SPEED_OF_LIGHT_M_S;
        finite(next, "travel time")?;
        if (next - dt).abs() <= tolerance_s {
            return Ok(Arrival {
                epoch: emission.shifted(next)?,
                travel_time_s: next,
                iterations,
            });
        }
        dt = next;
    }
    Err(Error::InvalidInput(
        "light-time iteration did not converge".into(),
    ))
}

#[derive(Debug, Clone, Copy, Default)]
pub struct DelayBudget {
    pub geometric_s: f64,
    pub medium_s: f64,
    pub clock_s: f64,
    pub gravitational_s: f64,
    pub instrumental_s: f64,
}
impl DelayBudget {
    pub fn total_s(self) -> Result<f64> {
        let terms = [
            self.geometric_s,
            self.medium_s,
            self.clock_s,
            self.gravitational_s,
            self.instrumental_s,
        ];
        for t in terms {
            finite(t, "delay term")?;
        }
        let total = terms.iter().sum();
        finite(total, "total delay")?;
        Ok(total)
    }
    /// Wrapped phase in radians, positive for positive arrival delay.
    pub fn phase_rad(self, frequency_hz: f64) -> Result<f64> {
        finite(frequency_hz, "frequency")?;
        if frequency_hz <= 0.0 {
            return Err(Error::InvalidInput("frequency must be positive".into()));
        }
        let cycles = self.total_s()? * frequency_hz;
        finite(cycles, "phase cycles")?;
        Ok(cycles.rem_euclid(1.0) * 2.0 * std::f64::consts::PI)
    }
}

/// Excess delay of a homogeneous, non-dispersive medium over vacuum.
pub fn medium_delay(path_length_m: f64, refractive_index: f64) -> Result<f64> {
    finite(path_length_m, "path length")?;
    finite(refractive_index, "refractive index")?;
    if path_length_m < 0.0 || refractive_index < 1.0 {
        return Err(Error::InvalidInput(
            "path must be nonnegative and refractive index >=1".into(),
        ));
    }
    let delay = path_length_m * (refractive_index - 1.0) / SPEED_OF_LIGHT_M_S;
    finite(delay, "medium delay")?;
    Ok(delay)
}

/// First-order weak-field Shapiro light-time term from a stationary point mass.
/// Endpoints and mass position share a common frame. Invalid at/through the mass.
/// Caller must separately reject occulted rays and supply body radii.
pub fn shapiro_delay(
    emitter: Vec3,
    receiver: Vec3,
    mass_position: Vec3,
    mu_m3_s2: f64,
) -> Result<f64> {
    emitter.validate()?;
    receiver.validate()?;
    mass_position.validate()?;
    finite(mu_m3_s2, "mu")?;
    if mu_m3_s2 < 0.0 {
        return Err(Error::InvalidInput("mu must be nonnegative".into()));
    }
    let r1 = (emitter - mass_position).norm();
    let r2 = (receiver - mass_position).norm();
    let r = (receiver - emitter).norm();
    let denominator = r1 + r2 - r;
    if r1 == 0.0 || r2 == 0.0 || denominator <= 0.0 {
        return Err(Error::InvalidInput("Shapiro ray is singular".into()));
    }
    let delay = 2.0 * mu_m3_s2 / SPEED_OF_LIGHT_M_S.powi(3) * (2.0 * r / denominator).ln_1p();
    finite(delay, "Shapiro delay")?;
    Ok(delay)
}

/// Clock reading minus coordinate time, in seconds; drift is s/s.
/// Calibration and time-scale conversion remain the caller's responsibility.
pub struct ClockModel {
    pub reference_epoch: Epoch,
    pub offset_s: f64,
    pub drift_s_per_s: f64,
    pub drift_rate_s_per_s2: f64,
}
impl ClockModel {
    pub fn offset_at(&self, epoch: Epoch) -> Result<f64> {
        for v in [self.offset_s, self.drift_s_per_s, self.drift_rate_s_per_s2] {
            finite(v, "clock coefficient")?;
        }
        let dt = epoch.duration_since(self.reference_epoch);
        let offset =
            self.offset_s + self.drift_s_per_s * dt + 0.5 * self.drift_rate_s_per_s2 * dt * dt;
        finite(offset, "clock offset")?;
        Ok(offset)
    }
}

/// Expected complex coherence for Gaussian delay uncertainty at one frequency.
/// Unit source visibility, positive phase convention matching DelayBudget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coherence {
    pub real: f64,
    pub imaginary: f64,
}
pub fn expected_coherence(
    mean_delay_s: f64,
    std_delay_s: f64,
    frequency_hz: f64,
) -> Result<Coherence> {
    finite(mean_delay_s, "mean delay")?;
    finite(std_delay_s, "delay uncertainty")?;
    finite(frequency_hz, "frequency")?;
    if std_delay_s < 0.0 || frequency_hz <= 0.0 {
        return Err(Error::InvalidInput(
            "require uncertainty>=0 and frequency>0".into(),
        ));
    }
    let cycles = frequency_hz * mean_delay_s;
    let std_phase = 2.0 * std::f64::consts::PI * frequency_hz * std_delay_s;
    finite(cycles, "phase cycles")?;
    finite(std_phase, "phase uncertainty")?;
    let phase = 2.0 * std::f64::consts::PI * cycles.rem_euclid(1.0);
    let amplitude = (-0.5 * std_phase * std_phase).exp();
    Ok(Coherence {
        real: amplitude * phase.cos(),
        imaginary: amplitude * phase.sin(),
    })
}

/// Special-relativistic frequency ratio observed/emitted, with velocities and
/// propagation direction in one inertial frame. Direction points along photon travel.
/// Does not include gravitational redshift, cosmological expansion, or moving media.
pub fn doppler_ratio(
    emitter_velocity_m_s: Vec3,
    receiver_velocity_m_s: Vec3,
    propagation_direction: Vec3,
) -> Result<f64> {
    emitter_velocity_m_s.validate()?;
    receiver_velocity_m_s.validate()?;
    let direction = propagation_direction.normalized()?;
    let be = emitter_velocity_m_s * (1.0 / SPEED_OF_LIGHT_M_S);
    let br = receiver_velocity_m_s * (1.0 / SPEED_OF_LIGHT_M_S);
    if be.norm() >= 1.0 || br.norm() >= 1.0 {
        return Err(Error::InvalidInput(
            "Doppler velocities must be subluminal".into(),
        ));
    }
    let ratio = ((1.0 - be.dot(be)) / (1.0 - br.dot(br))).sqrt() * (1.0 - direction.dot(br))
        / (1.0 - direction.dot(be));
    finite(ratio, "Doppler ratio")?;
    Ok(ratio)
}

/// A phase-centre ray can be excluded by an opaque spherical body. Coordinates
/// share a common frame; the closed line segment includes its endpoints.
pub fn ray_intersects_sphere(
    emitter: Vec3,
    receiver: Vec3,
    centre: Vec3,
    radius_m: f64,
) -> Result<bool> {
    emitter.validate()?;
    receiver.validate()?;
    centre.validate()?;
    finite(radius_m, "radius")?;
    if radius_m <= 0.0 {
        return Err(Error::InvalidInput("sphere radius must be positive".into()));
    }
    let ray = receiver - emitter;
    let length = ray.norm();
    if length == 0.0 {
        return Ok((emitter - centre).norm() <= radius_m);
    }
    let direction = ray * (1.0 / length);
    let along = (centre - emitter).dot(direction).clamp(0.0, length);
    Ok((emitter + direction * along - centre).norm() <= radius_m)
}
