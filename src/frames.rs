//! Epoch-aware frame tree for arbitrary bodies and star systems.
//!
//! The root is a user-defined common inertial reference, not automatically ICRF.
//! Transformations are Newtonian rigid transforms evaluated at one coordinate time.
use crate::coordinates::Ellipsoid;
use crate::{Error, Result, finite};
use std::ops::{Add, Mul, Neg, Sub};
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Vec3 {
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self> {
        finite(x, "x")?;
        finite(y, "y")?;
        finite(z, "z")?;
        Ok(Self { x, y, z })
    }
    pub fn dot(self, b: Self) -> f64 {
        self.x * b.x + self.y * b.y + self.z * b.z
    }
    pub fn cross(self, b: Self) -> Self {
        Self {
            x: self.y * b.z - self.z * b.y,
            y: self.z * b.x - self.x * b.z,
            z: self.x * b.y - self.y * b.x,
        }
    }
    pub fn norm(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }
    pub fn normalized(self) -> Result<Self> {
        let n = self.norm();
        if !n.is_finite() || n == 0.0 {
            return Err(Error::InvalidInput(
                "direction must be finite and nonzero".into(),
            ));
        }
        Ok(self * (1.0 / n))
    }
    pub fn validate(self) -> Result<()> {
        finite(self.x, "x")?;
        finite(self.y, "y")?;
        finite(self.z, "z")
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            x: self.x + b.x,
            y: self.y + b.y,
            z: self.z + b.z,
        }
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            x: self.x - b.x,
            y: self.y - b.y,
            z: self.z - b.z,
        }
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self {
            x: self.x * k,
            y: self.y * k,
            z: self.z * k,
        }
    }
}
impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}

/// Continuous coordinate seconds relative to a caller-chosen origin.
/// A single time scale must be used throughout a graph (e.g. TDB for ephemerides).
/// UTC/leap-second and relativistic time-scale conversions are intentionally external.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Epoch {
    whole: i64,
    fraction: f64,
}
impl Epoch {
    pub fn new(seconds: f64) -> Result<Self> {
        Self::from_parts(0, seconds)
    }
    pub fn from_parts(whole_seconds: i64, fractional_seconds: f64) -> Result<Self> {
        finite(fractional_seconds, "epoch")?;
        let shift = fractional_seconds.floor();
        if !(-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&shift) {
            return Err(Error::InvalidInput("epoch overflow".into()));
        }
        let whole = whole_seconds
            .checked_add(shift as i64)
            .ok_or_else(|| Error::InvalidInput("epoch overflow".into()))?;
        let fraction = fractional_seconds - shift;
        if fraction >= 1.0 {
            return Ok(Self {
                whole: whole
                    .checked_add(1)
                    .ok_or_else(|| Error::InvalidInput("epoch overflow".into()))?,
                fraction: 0.0,
            });
        }
        Ok(Self { whole, fraction })
    }
    /// Approximate scalar seconds, for display/interchange only.
    /// Use duration_since for calculations, and from_parts for fine epoch input.
    pub fn seconds(self) -> f64 {
        self.whole as f64 + self.fraction
    }
    pub fn whole_seconds(self) -> i64 {
        self.whole
    }
    pub fn fractional_seconds(self) -> f64 {
        self.fraction
    }
    pub fn duration_since(self, other: Self) -> f64 {
        (i128::from(self.whole) - i128::from(other.whole)) as f64 + (self.fraction - other.fraction)
    }
    pub fn shifted(self, seconds: f64) -> Result<Self> {
        finite(seconds, "epoch shift")?;
        let shift = Self::new(seconds)?;
        let whole = self
            .whole
            .checked_add(shift.whole)
            .ok_or_else(|| Error::InvalidInput("epoch overflow".into()))?;
        Self::from_parts(whole, self.fraction + shift.fraction)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    pub position_m: Vec3,
    pub velocity_m_s: Vec3,
}
impl State {
    pub fn stationary(position_m: Vec3) -> Self {
        Self {
            position_m,
            velocity_m_s: Vec3::ZERO,
        }
    }
}

/// Proper orthonormal rotation from local axes to parent axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotation {
    rows: [[f64; 3]; 3],
}
impl Rotation {
    pub const IDENTITY: Self = Self {
        rows: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };
    pub fn new(rows: [[f64; 3]; 3]) -> Result<Self> {
        let vectors = rows.map(|r| Vec3 {
            x: r[0],
            y: r[1],
            z: r[2],
        });
        for (i, a) in vectors.iter().enumerate() {
            a.validate()?;
            for (j, b) in vectors.iter().enumerate() {
                if (a.dot(*b) - if i == j { 1.0 } else { 0.0 }).abs() > 1e-10 {
                    return Err(Error::InvalidInput("rotation must be orthonormal".into()));
                }
            }
        }
        if (vectors[0].cross(vectors[1]).dot(vectors[2]) - 1.0).abs() > 1e-10 {
            return Err(Error::InvalidInput(
                "rotation must preserve handedness".into(),
            ));
        }
        Ok(Self { rows })
    }
    pub fn axis_angle(axis: Vec3, angle_rad: f64) -> Result<Self> {
        finite(angle_rad, "rotation angle")?;
        let a = axis.normalized()?;
        let (s, c) = angle_rad.sin_cos();
        let t = 1.0 - c;
        Self::new([
            [
                c + a.x * a.x * t,
                a.x * a.y * t - a.z * s,
                a.x * a.z * t + a.y * s,
            ],
            [
                a.y * a.x * t + a.z * s,
                c + a.y * a.y * t,
                a.y * a.z * t - a.x * s,
            ],
            [
                a.z * a.x * t - a.y * s,
                a.z * a.y * t + a.x * s,
                c + a.z * a.z * t,
            ],
        ])
    }
    pub fn apply(self, v: Vec3) -> Vec3 {
        let r = self.rows;
        Vec3 {
            x: r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z,
            y: r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
            z: r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z,
        }
    }
    pub fn inverse(self) -> Self {
        let r = self.rows;
        Self {
            rows: [
                [r[0][0], r[1][0], r[2][0]],
                [r[0][1], r[1][1], r[2][1]],
                [r[0][2], r[1][2], r[2][2]],
            ],
        }
    }
    pub fn compose(self, b: Self) -> Self {
        let mut r = [[0.0; 3]; 3];
        for (i, row) in r.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = (0..3).map(|k| self.rows[i][k] * b.rows[k][j]).sum();
            }
        }
        Self { rows: r }
    }
}

/// Position, velocity, orientation, and angular velocity of a frame in its parent.
#[derive(Debug, Clone, Copy)]
pub struct Pose {
    pub origin: State,
    pub rotation: Rotation,
    pub angular_velocity_rad_s: Vec3,
}
impl Pose {
    pub const IDENTITY: Self = Self {
        origin: State {
            position_m: Vec3::ZERO,
            velocity_m_s: Vec3::ZERO,
        },
        rotation: Rotation::IDENTITY,
        angular_velocity_rad_s: Vec3::ZERO,
    };
    pub fn apply(self, state: State) -> State {
        let offset = self.rotation.apply(state.position_m);
        State {
            position_m: self.origin.position_m + offset,
            velocity_m_s: self.origin.velocity_m_s
                + self.rotation.apply(state.velocity_m_s)
                + self.angular_velocity_rad_s.cross(offset),
        }
    }
    pub fn unapply(self, state: State) -> State {
        let offset = state.position_m - self.origin.position_m;
        State {
            position_m: self.rotation.inverse().apply(offset),
            velocity_m_s: self.rotation.inverse().apply(
                state.velocity_m_s
                    - self.origin.velocity_m_s
                    - self.angular_velocity_rad_s.cross(offset),
            ),
        }
    }
    fn compose(self, child: Self) -> Self {
        Self {
            origin: self.apply(child.origin),
            rotation: self.rotation.compose(child.rotation),
            angular_velocity_rad_s: self.angular_velocity_rad_s
                + self.rotation.apply(child.angular_velocity_rad_s),
        }
    }
    fn validate(self) -> Result<Self> {
        self.origin.position_m.validate()?;
        self.origin.velocity_m_s.validate()?;
        self.angular_velocity_rad_s.validate()?;
        Ok(self)
    }
}

/// Implement this for authoritative ephemerides and body-orientation models.
pub trait PoseProvider: Send + Sync {
    fn pose(&self, epoch: Epoch) -> Result<Pose>;
}
impl PoseProvider for Pose {
    fn pose(&self, _: Epoch) -> Result<Pose> {
        self.validate()
    }
}

pub struct LinearMotion {
    pub reference_epoch: Epoch,
    pub origin: State,
    pub initial_rotation: Rotation,
    pub angular_velocity_rad_s: Vec3,
}
impl PoseProvider for LinearMotion {
    fn pose(&self, epoch: Epoch) -> Result<Pose> {
        let dt = epoch.duration_since(self.reference_epoch);
        let speed = self.angular_velocity_rad_s.norm();
        let rotation = if speed == 0.0 {
            self.initial_rotation
        } else {
            Rotation::axis_angle(self.angular_velocity_rad_s, speed * dt)?
                .compose(self.initial_rotation)
        };
        Pose {
            origin: State {
                position_m: self.origin.position_m + self.origin.velocity_m_s * dt,
                velocity_m_s: self.origin.velocity_m_s,
            },
            rotation,
            angular_velocity_rad_s: self.angular_velocity_rad_s,
        }
        .validate()
    }
}

/// Two-body bound Kepler orbit about the parent origin. Angles are radians.
/// No perturbations, light-time, precession, or relativistic corrections are implicit.
pub struct KeplerOrbit {
    pub reference_epoch: Epoch,
    pub semi_major_axis_m: f64,
    pub eccentricity: f64,
    pub gravitational_parameter_m3_s2: f64,
    pub mean_anomaly_rad: f64,
    pub orbital_plane: Rotation,
}
impl PoseProvider for KeplerOrbit {
    fn pose(&self, epoch: Epoch) -> Result<Pose> {
        for (v, n) in [
            (self.semi_major_axis_m, "semi-major axis"),
            (self.eccentricity, "eccentricity"),
            (self.gravitational_parameter_m3_s2, "mu"),
            (self.mean_anomaly_rad, "mean anomaly"),
        ] {
            finite(v, n)?;
        }
        if self.semi_major_axis_m <= 0.0
            || self.gravitational_parameter_m3_s2 <= 0.0
            || !(0.0..1.0).contains(&self.eccentricity)
        {
            return Err(Error::InvalidInput(
                "Kepler model requires a>0, mu>0, and 0<=e<1".into(),
            ));
        }
        let a = self.semi_major_axis_m;
        let e = self.eccentricity;
        let n = (self.gravitational_parameter_m3_s2 / a.powi(3)).sqrt();
        let m = (self.mean_anomaly_rad + n * epoch.duration_since(self.reference_epoch))
            .rem_euclid(2.0 * std::f64::consts::PI);
        // Monotonic bisection avoids the high-eccentricity Newton divergence.
        let (mut lo, mut hi) = (0.0, 2.0 * std::f64::consts::PI);
        for _ in 0..60 {
            let mid = (lo + hi) / 2.0;
            if mid - e * mid.sin() < m {
                lo = mid
            } else {
                hi = mid
            }
        }
        let anomaly = (lo + hi) / 2.0;
        let (s, c) = anomaly.sin_cos();
        let q = (1.0 - e * e).sqrt();
        let rate = n / (1.0 - e * c);
        Pose {
            origin: State {
                position_m: self.orbital_plane.apply(Vec3 {
                    x: a * (c - e),
                    y: a * q * s,
                    z: 0.0,
                }),
                velocity_m_s: self.orbital_plane.apply(Vec3 {
                    x: -a * s * rate,
                    y: a * q * c * rate,
                    z: 0.0,
                }),
            },
            rotation: Rotation::IDENTITY,
            angular_velocity_rad_s: Vec3::ZERO,
        }
        .validate()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FrameId {
    graph: u64,
    index: usize,
}
struct Node {
    name: String,
    parent: Option<FrameId>,
    provider: Box<dyn PoseProvider>,
}
pub struct FrameGraph {
    id: u64,
    nodes: Vec<Node>,
    time_reference: Option<crate::time::TimeReference>,
}
impl FrameGraph {
    pub fn new(root_name: impl Into<String>) -> Self {
        static NEXT_ID: AtomicU64 = AtomicU64::new(1);
        Self {
            id: NEXT_ID.fetch_add(1, Ordering::Relaxed),
            nodes: vec![Node {
                name: root_name.into(),
                parent: None,
                provider: Box::new(Pose::IDENTITY),
            }],
            time_reference: None,
        }
    }
    pub fn with_time_reference(
        root_name: impl Into<String>,
        reference: crate::time::TimeReference,
    ) -> Self {
        let mut graph = Self::new(root_name);
        graph.time_reference = Some(reference);
        graph
    }
    pub fn time_reference(&self) -> Option<crate::time::TimeReference> {
        self.time_reference
    }
    pub fn epoch_at(&self, timestamp: crate::time::Timestamp) -> Result<Epoch> {
        self.time_reference
            .ok_or_else(|| Error::InvalidInput("graph has no declared time reference".into()))?
            .epoch_at(timestamp)
    }
    pub fn root(&self) -> FrameId {
        FrameId {
            graph: self.id,
            index: 0,
        }
    }
    pub fn add(
        &mut self,
        name: impl Into<String>,
        parent: FrameId,
        provider: impl PoseProvider + 'static,
    ) -> Result<FrameId> {
        self.node(parent)?;
        let name = name.into();
        if name.trim().is_empty() || self.nodes.iter().any(|n| n.name == name) {
            return Err(Error::InvalidInput(
                "frame names must be nonempty and unique".into(),
            ));
        }
        let id = FrameId {
            graph: self.id,
            index: self.nodes.len(),
        };
        self.nodes.push(Node {
            name,
            parent: Some(parent),
            provider: Box::new(provider),
        });
        Ok(id)
    }
    fn node(&self, id: FrameId) -> Result<&Node> {
        if id.graph != self.id {
            return Err(Error::InvalidInput(
                "frame belongs to a different graph".into(),
            ));
        }
        self.nodes
            .get(id.index)
            .ok_or_else(|| Error::InvalidInput("unknown frame".into()))
    }
    pub fn name(&self, id: FrameId) -> Result<&str> {
        Ok(&self.node(id)?.name)
    }
    pub fn pose_in_root(&self, frame: FrameId, epoch: Epoch) -> Result<Pose> {
        let mut pose = Pose::IDENTITY;
        let mut current = frame;
        loop {
            let node = self.node(current)?;
            pose = node.provider.pose(epoch)?.validate()?.compose(pose);
            match node.parent {
                Some(parent) => current = parent,
                None => return pose.validate(),
            }
        }
    }
    pub fn transform(
        &self,
        state: State,
        from: FrameId,
        to: FrameId,
        epoch: Epoch,
    ) -> Result<State> {
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        // Traverse only to the lowest common ancestor. Shared astronomical
        // translations cancel symbolically, preserving small local baselines.
        let mut from_path = vec![from];
        let mut current = from;
        while let Some(parent) = self.node(current)?.parent {
            from_path.push(parent);
            current = parent;
        }
        let mut to_path = vec![];
        current = to;
        while !from_path.contains(&current) {
            let node = self.node(current)?;
            to_path.push(current);
            current = node
                .parent
                .ok_or_else(|| Error::InvalidInput("frames have no common ancestor".into()))?;
        }
        self.node(current)?;
        let common = current;
        let mut result = state;
        current = from;
        while current != common {
            let node = self.node(current)?;
            result = node.provider.pose(epoch)?.validate()?.apply(result);
            current = node.parent.expect("validated ancestor path");
        }
        for id in to_path.iter().rev() {
            result = self
                .node(*id)?
                .provider
                .pose(epoch)?
                .validate()?
                .unapply(result);
        }
        result.position_m.validate()?;
        result.velocity_m_s.validate()?;
        Ok(result)
    }
}

/// Associates an ellipsoid and gravitational parameter with a body-fixed frame.
/// Arbitrarily many planets and stars may coexist beneath a common root.
pub struct Body {
    pub name: String,
    pub fixed_frame: FrameId,
    pub ellipsoid: Ellipsoid,
    pub gravitational_parameter_m3_s2: f64,
}
impl Body {
    pub fn surface_state(&self, point: crate::coordinates::Geodetic) -> State {
        let p = self.ellipsoid.to_cartesian(point);
        State::stationary(Vec3 {
            x: p.x,
            y: p.y,
            z: p.z,
        })
    }
}
