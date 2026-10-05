//! Geospatial primitives and explicit, source-backed physical expectation models.
//!
//! Angles at API boundaries are degrees; lengths are metres; RF powers are watts
//! unless explicitly named `dbm`. All coordinate heights are ellipsoidal unless
//! a vertical conversion is explicitly requested. See the README for model limits.
//!
//! ```
//! use geoid_atlas::coordinates::{Ellipsoid, Geodetic};
//! use geoid_atlas::frames::{Epoch, FrameGraph, Pose, State, Vec3};
//! use geoid_atlas::interferometry::far_field;
//!
//! let mut frames = FrameGraph::new("common inertial reference");
//! let root = frames.root();
//! let body = frames.add("body-fixed", root, Pose::IDENTITY)?;
//! let site = Ellipsoid::WGS84.to_cartesian(Geodetic::new(0.0, 0.0, -100.0)?);
//! let a = frames.transform(State::stationary(Vec3::new(site.x, site.y, site.z)?),
//!                          body, root, Epoch::new(0.0)?)?;
//! let b = State::stationary(a.position_m + Vec3::new(100.0, 0.0, 0.0)?);
//! let delay = far_field(a, b, Vec3::new(1.0, 0.0, 0.0)?)?;
//! assert!(delay.seconds < 0.0); // receiver 2 is closer to the source
//! # Ok::<(), geoid_atlas::Error>(())
//! ```

pub mod atlas;
#[cfg(any(feature = "proj-backend", feature = "gdal-backend"))]
pub mod backends;
pub mod calibration;
pub mod coordinates;
pub mod echo;
pub mod environment;
pub mod estimation;
pub mod formats;
pub mod frames;
pub mod fusion;
pub mod geodesics;
pub mod gravity;
pub mod interferometry;
pub mod maidenhead;
pub mod observability;
pub mod raster;
pub mod reference;
pub mod relativity;
pub mod rf;
pub mod spacetime;
pub mod time;
pub mod units;

use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Backend { backend: String, message: String },
    InvalidInput(String),
    OutsideCoverage,
    NoData,
    Parse(String),
    DuplicateSource(String),
    UnknownSource(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Backend { backend, message } => write!(f, "{backend} backend: {message}"),
            Self::InvalidInput(s) => write!(f, "invalid input: {s}"),
            Self::OutsideCoverage => f.write_str("point lies outside coverage"),
            Self::NoData => f.write_str("no data at point"),
            Self::Parse(s) => write!(f, "parse error: {s}"),
            Self::DuplicateSource(s) => write!(f, "duplicate source: {s}"),
            Self::UnknownSource(s) => write!(f, "unknown source: {s}"),
        }
    }
}
impl std::error::Error for Error {}
pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!("{name} must be finite")))
    }
}
