//! Karney surface geodesics on an oblate reference ellipsoid.
//!
//! Coordinates use geodetic latitude. Heights are ignored: these are surface
//! distances, not chord lengths, terrain-following paths or sky-wave ray lengths.
//! Azimuths are degrees clockwise from north; the endpoint azimuth is **forward**
//! along the path, not the bearing back to its origin. Longitudes are normalized
//! to [-180, 180). At coincident or antipodal points the azimuth/path may not be
//! unique; the solver selects GeographicLib's canonical solution.
//!
//! The pure Rust `geographiclib-rs` backend implements Karney's sixth-order
//! algorithm, including the near-antipodal inverse solve. WGS84 results have
//! round-off-level accuracy; accuracy decreases for larger flattening. This API
//! accepts only flattening <= 0.02, rather than silently using the series on
//! highly flattened bodies. No native GIS library or geoid grid is required.
//!
//! ```
//! use geoid_atlas::coordinates::Geodetic;
//! use geoid_atlas::geodesics::{SurfaceGeodesic, ground_distance_m};
//! let london = Geodetic::new(51.5, -0.12, 0.0)?;
//! let wellington = Geodetic::new(-41.3, 174.78, 0.0)?;
//! let path = SurfaceGeodesic::wgs84().path(london, wellington);
//! let midpoint = path.midpoint();
//! assert!((ground_distance_m(london, midpoint) - path.inverse().distance_m / 2.0).abs() < 1e-7);
//! # Ok::<(), geoid_atlas::Error>(())
//! ```

use crate::coordinates::{Ellipsoid, Geodetic};
use crate::{Error, Result, finite};
use geographiclib_rs::{DirectGeodesic, Geodesic, InverseGeodesic};

/// IUGG mean Earth radius, matching the conventional FT8 great-circle model.
pub const MEAN_EARTH_RADIUS_M: f64 = 6_371_008.8;

/// Shortest surface path between two points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodesicInverse {
    pub distance_m: f64,
    pub initial_azimuth_deg: f64,
    pub final_azimuth_deg: f64,
    /// Angle on the auxiliary sphere; this is not distance / mean Earth radius.
    pub auxiliary_arc_deg: f64,
}

/// Surface endpoint and forward heading after travelling a signed distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodesicDirect {
    /// Height is zero relative to the chosen reference ellipsoid.
    pub point: Geodetic,
    pub final_azimuth_deg: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SurfaceGeodesic {
    ellipsoid: Ellipsoid,
    solver: Geodesic,
}

impl SurfaceGeodesic {
    /// Cached WGS84 solver; no fallible construction is needed.
    pub fn wgs84() -> Self {
        Self {
            ellipsoid: Ellipsoid::WGS84,
            solver: Geodesic::wgs84(),
        }
    }

    /// Configure a sphere or oblate ellipsoid with flattening in [0, 0.02].
    /// The squared radius must be a normal finite f64 to avoid numerical overflow
    /// or underflow in the backend. This easily includes physical planetary sizes.
    pub fn new(ellipsoid: Ellipsoid) -> Result<Self> {
        let radius_squared = ellipsoid.equatorial_radius_m().powi(2);
        if ellipsoid.flattening() > 0.02 || !radius_squared.is_normal() {
            return Err(Error::InvalidInput(
                "surface geodesics require flattening <= 0.02 and a normal finite squared radius"
                    .into(),
            ));
        }
        Ok(Self {
            ellipsoid,
            solver: Geodesic::new(ellipsoid.equatorial_radius_m(), ellipsoid.flattening()),
        })
    }

    pub fn ellipsoid(self) -> Ellipsoid {
        self.ellipsoid
    }

    /// Solve the shortest inverse path; input heights are ignored.
    pub fn inverse(self, start: Geodetic, end: Geodetic) -> GeodesicInverse {
        let (distance_m, initial_azimuth_deg, final_azimuth_deg, auxiliary_arc_deg) =
            self.solver.inverse(
                start.latitude_deg(),
                start.longitude_deg(),
                end.latitude_deg(),
                end.longitude_deg(),
            );
        GeodesicInverse {
            distance_m,
            initial_azimuth_deg,
            final_azimuth_deg,
            auxiliary_arc_deg,
        }
    }

    /// Follow an azimuth by a signed distance; negative distances travel backward.
    /// Azimuth and distance must be finite. The endpoint has ellipsoidal height 0.
    /// Very large distances that overflow the solver are rejected.
    pub fn direct(
        self,
        start: Geodetic,
        azimuth_deg: f64,
        distance_m: f64,
    ) -> Result<GeodesicDirect> {
        finite(azimuth_deg, "azimuth")?;
        finite(distance_m, "surface distance")?;
        let (latitude, longitude, final_azimuth_deg) = self.solver.direct(
            start.latitude_deg(),
            start.longitude_deg(),
            azimuth_deg,
            distance_m,
        );
        finite(final_azimuth_deg, "computed azimuth")?;
        Ok(GeodesicDirect {
            point: Geodetic::new(latitude, longitude, 0.0)?,
            final_azimuth_deg,
        })
    }

    /// The canonical shortest geodesic, with interpolation by surface distance.
    pub fn path(self, start: Geodetic, end: Geodetic) -> GeodesicPath {
        GeodesicPath {
            geodesic: self,
            start,
            inverse: self.inverse(start, end),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GeodesicPath {
    geodesic: SurfaceGeodesic,
    start: Geodetic,
    inverse: GeodesicInverse,
}

impl GeodesicPath {
    pub fn inverse(self) -> GeodesicInverse {
        self.inverse
    }

    /// Point and forward heading at distance in [0, path length].
    /// Use `SurfaceGeodesic::direct` for extrapolation beyond the endpoints.
    pub fn position(self, distance_m: f64) -> Result<GeodesicDirect> {
        finite(distance_m, "path distance")?;
        if !(0.0..=self.inverse.distance_m).contains(&distance_m) {
            return Err(Error::InvalidInput(
                "distance lies outside surface path".into(),
            ));
        }
        self.geodesic
            .direct(self.start, self.inverse.initial_azimuth_deg, distance_m)
    }

    /// Fraction in [0, 1] of the total surface distance; output height is zero.
    pub fn point_at_fraction(self, fraction: f64) -> Result<Geodetic> {
        finite(fraction, "path fraction")?;
        if !(0.0..=1.0).contains(&fraction) {
            return Err(Error::InvalidInput("path fraction must be in [0,1]".into()));
        }
        Ok(self.position(fraction * self.inverse.distance_m)?.point)
    }

    /// Halfway along the surface geodesic, distinct from a Cartesian chord midpoint.
    pub fn midpoint(self) -> Geodetic {
        self.point_at_fraction(0.5)
            .expect("validated geodesic has a finite midpoint")
    }
}

/// WGS84 Karney distance in metres, ignoring heights. This deliberately replaces
/// Lambert's approximation; existing calibration values can therefore shift.
pub fn ground_distance_m(start: Geodetic, end: Geodetic) -> f64 {
    SurfaceGeodesic::wgs84().inverse(start, end).distance_m
}

/// Spherical great-circle distance using R = 6,371,008.8 m.
/// Geodetic latitudes are used as spherical latitudes for compatibility with
/// conventional FT8 ground-path models; heights are ignored.
pub fn great_circle_m(start: Geodetic, end: Geodetic) -> f64 {
    let (s1, c1) = start.latitude_deg().to_radians().sin_cos();
    let (s2, c2) = end.latitude_deg().to_radians().sin_cos();
    let delta = (end.longitude_deg() - start.longitude_deg()).to_radians();
    let (sd, cd) = delta.sin_cos();
    // atan2(cross, dot) stays stable at both small and near-antipodal angles.
    let cross = (c2 * sd).hypot(c1 * s2 - s1 * c2 * cd);
    let dot = s1 * s2 + c1 * c2 * cd;
    MEAN_EARTH_RADIUS_M * cross.atan2(dot)
}
