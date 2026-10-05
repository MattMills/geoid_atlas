//! Normal Earth gravity, geoid vertical conversions, and multi-body point-mass fields.
use crate::coordinates::{Geodetic, WGS84_A, WGS84_F};
use crate::frames::Vec3;
use crate::raster::Raster;
use crate::{Error, Result, finite};

/// Somigliana normal gravity on the WGS84 reference ellipsoid, m/s².
/// This is apparent normal gravity (including rotation), not a measured anomaly.
pub fn wgs84_normal_gravity(latitude_deg: f64) -> Result<f64> {
    finite(latitude_deg, "latitude")?;
    if !(-90.0..=90.0).contains(&latitude_deg) {
        return Err(Error::InvalidInput("latitude outside [-90,90]".into()));
    }
    let sin2 = latitude_deg.to_radians().sin().powi(2);
    let e2 = WGS84_F * (2.0 - WGS84_F);
    Ok(9.780_325_335_9 * (1.0 + 0.001_931_852_652_41 * sin2) / (1.0 - e2 * sin2).sqrt())
}

/// Low-altitude first-order free-air approximation; valid only near Earth's surface.
pub fn wgs84_free_air_gravity(point: Geodetic) -> Result<f64> {
    if point.height_m().abs() > 10_000.0 {
        return Err(Error::InvalidInput(
            "free-air approximation limited to +/-10 km".into(),
        ));
    }
    Ok(wgs84_normal_gravity(point.latitude_deg())? - 3.086e-6 * point.height_m())
}

pub trait GeoidModel {
    fn undulation_m(&self, point: Geodetic) -> Result<f64>;
}
impl GeoidModel for Raster<f64> {
    fn undulation_m(&self, point: Geodetic) -> Result<f64> {
        self.bilinear(point)
    }
}

/// h = H + N. A geoid model must match the body's horizontal/vertical reference.
pub fn ellipsoidal_height(orthometric_height_m: f64, undulation_m: f64) -> Result<f64> {
    finite(orthometric_height_m, "orthometric height")?;
    finite(undulation_m, "undulation")?;
    let h = orthometric_height_m + undulation_m;
    finite(h, "ellipsoidal height")?;
    Ok(h)
}
pub fn orthometric_height(ellipsoidal_height_m: f64, undulation_m: f64) -> Result<f64> {
    finite(ellipsoidal_height_m, "ellipsoidal height")?;
    finite(undulation_m, "undulation")?;
    let h = ellipsoidal_height_m - undulation_m;
    finite(h, "orthometric height")?;
    Ok(h)
}

#[derive(Debug, Clone, Copy)]
pub struct PointMass {
    pub position_m: Vec3,
    pub gravitational_parameter_m3_s2: f64,
}
#[derive(Debug, Clone, Copy, Default)]
pub struct GravityField {
    pub potential_m2_s2: f64,
    pub acceleration_m_s2: Vec3,
}

/// Newtonian superposition in a common frame at a common epoch.
/// Potential is -mu/r, zero at infinity. No centrifugal term is included.
pub fn point_mass_field(point: Vec3, masses: &[PointMass]) -> Result<GravityField> {
    point.validate()?;
    let mut result = GravityField::default();
    for mass in masses {
        mass.position_m.validate()?;
        let mu = mass.gravitational_parameter_m3_s2;
        finite(mu, "mu")?;
        if mu < 0.0 {
            return Err(Error::InvalidInput("mu must be nonnegative".into()));
        }
        if mu == 0.0 {
            continue;
        }
        let delta = mass.position_m - point;
        let r = delta.norm();
        if r == 0.0 {
            return Err(Error::InvalidInput(
                "point-mass field is singular at the mass".into(),
            ));
        }
        result.potential_m2_s2 -= mu / r;
        result.acceleration_m_s2 = result.acceleration_m_s2 + delta * (mu / r / r / r);
    }
    finite(result.potential_m2_s2, "potential")?;
    result.acceleration_m_s2.validate()?;
    Ok(result)
}

/// WGS84 equatorial reference radius, provided for callers supplying gravity coefficients.
pub const EARTH_REFERENCE_RADIUS_M: f64 = WGS84_A;
