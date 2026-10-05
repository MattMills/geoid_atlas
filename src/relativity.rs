//! Special-relativistic event/velocity transforms between inertial observers.
//! Kept separate from the Newtonian FrameGraph: a boost changes simultaneity.
//! No curved metric, barycentric time-scale conversion, gravitational light bending,
//! or proper-time clock synchronization is implied by these local flat-space tools.
use crate::frames::{Epoch, State, Vec3};
use crate::interferometry::SPEED_OF_LIGHT_M_S as C;
use crate::{Error, Result, finite};

#[derive(Debug, Clone, Copy)]
pub struct Event {
    pub epoch: Epoch,
    pub position_m: Vec3,
}
/// SI electromagnetic fields at one event. A boost mixes E and B; a magnetic
/// field alone cannot be transformed relativistically without an electric model.
#[derive(Debug, Clone, Copy)]
pub struct ElectromagneticField {
    pub electric_v_m: Vec3,
    pub magnetic_t: Vec3,
}
/// Standard boost with parallel axes: target origin moves at +velocity in source.
/// The supplied source origin event maps to target position/time zero. Use nearby
/// origins and split epochs; subtraction cannot restore already-rounded coordinates.
#[derive(Debug, Clone, Copy)]
pub struct LorentzBoost {
    source_origin: Event,
    velocity: Vec3,
    direction: Vec3,
    gamma: f64,
    gamma_minus_one: f64,
}
impl LorentzBoost {
    /// Transform fields to the same observer as apply(event), in parallel axes.
    pub fn electromagnetic_field(
        self,
        field: ElectromagneticField,
    ) -> Result<ElectromagneticField> {
        field.electric_v_m.validate()?;
        field.magnetic_t.validate()?;
        let electric_v_m = (field.electric_v_m + self.velocity.cross(field.magnetic_t))
            * self.gamma
            - self.direction * (self.gamma_minus_one * self.direction.dot(field.electric_v_m));
        let magnetic_t = (field.magnetic_t
            - self.velocity.cross(field.electric_v_m) * (1.0 / C / C))
            * self.gamma
            - self.direction * (self.gamma_minus_one * self.direction.dot(field.magnetic_t));
        electric_v_m.validate()?;
        magnetic_t.validate()?;
        Ok(ElectromagneticField {
            electric_v_m,
            magnetic_t,
        })
    }
    pub fn new(source_origin: Event, velocity_m_s: Vec3) -> Result<Self> {
        source_origin.position_m.validate()?;
        velocity_m_s.validate()?;
        let speed = velocity_m_s.norm();
        if speed >= C {
            return Err(Error::InvalidInput("boost must be subluminal".into()));
        }
        let beta = speed / C;
        let gamma = 1.0 / ((1.0 - beta) * (1.0 + beta)).sqrt();
        let gamma_minus_one = gamma * gamma * beta * beta / (gamma + 1.0);
        let direction = if speed == 0.0 {
            Vec3::ZERO
        } else {
            velocity_m_s * (1.0 / speed)
        };
        Ok(Self {
            source_origin,
            velocity: velocity_m_s,
            direction,
            gamma,
            gamma_minus_one,
        })
    }
    pub fn apply(self, event: Event) -> Result<Event> {
        event.position_m.validate()?;
        let x = event.position_m - self.source_origin.position_m;
        let dt = event.epoch.duration_since(self.source_origin.epoch);
        let seconds = self.gamma * (dt - self.velocity.dot(x) / C / C);
        let position_m = x + self.direction * (self.gamma_minus_one * self.direction.dot(x))
            - self.velocity * (self.gamma * dt);
        position_m.validate()?;
        Ok(Event {
            epoch: Epoch::new(seconds)?,
            position_m,
        })
    }
    pub fn unapply(self, event: Event) -> Result<Event> {
        event.position_m.validate()?;
        let dt = event.epoch.duration_since(Epoch::new(0.0)?);
        let x = event.position_m;
        let seconds = self.gamma * (dt + self.velocity.dot(x) / C / C);
        let position_m = self.source_origin.position_m
            + x
            + self.direction * (self.gamma_minus_one * self.direction.dot(x))
            + self.velocity * (self.gamma * dt);
        position_m.validate()?;
        Ok(Event {
            epoch: self.source_origin.epoch.shifted(seconds)?,
            position_m,
        })
    }
    /// Transform the state at its event, not a constant-time slice of an array.
    pub fn apply_state(self, epoch: Epoch, state: State) -> Result<(Epoch, State)> {
        state.velocity_m_s.validate()?;
        if state.velocity_m_s.norm() >= C {
            return Err(Error::InvalidInput("state must be subluminal".into()));
        }
        let event = self.apply(Event {
            epoch,
            position_m: state.position_m,
        })?;
        let denominator = self.gamma * (1.0 - self.velocity.dot(state.velocity_m_s) / C / C);
        finite(denominator, "velocity transform denominator")?;
        if denominator <= 0.0 {
            return Err(Error::InvalidInput("invalid velocity transform".into()));
        }
        let velocity_m_s = (state.velocity_m_s
            + self.direction * (self.gamma_minus_one * self.direction.dot(state.velocity_m_s))
            - self.velocity * self.gamma)
            * (1.0 / denominator);
        velocity_m_s.validate()?;
        Ok((
            event.epoch,
            State {
                position_m: event.position_m,
                velocity_m_s,
            },
        ))
    }
}
