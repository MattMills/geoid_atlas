//! WGS84 geodetic, Earth-centred Earth-fixed, local ENU, Web Mercator, and UTM.
use crate::frames::Vec3;
use crate::{Error, Result, finite};
use std::ops::{Add, Div, Mul, Sub};

pub const WGS84_A: f64 = 6_378_137.0;
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
const E2: f64 = WGS84_F * (2.0 - WGS84_F);

/// Oblate reference ellipsoid. Heights may be negative (subsurface).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    a: f64,
    flattening: f64,
}
impl Ellipsoid {
    pub const WGS84: Self = Self {
        a: WGS84_A,
        flattening: WGS84_F,
    };
    pub fn new(equatorial_radius_m: f64, flattening: f64) -> Result<Self> {
        finite(equatorial_radius_m, "equatorial radius")?;
        finite(flattening, "flattening")?;
        if equatorial_radius_m <= 0.0 || !(0.0..1.0).contains(&flattening) {
            return Err(Error::InvalidInput(
                "ellipsoid requires radius > 0 and flattening in [0,1)".into(),
            ));
        }
        Ok(Self {
            a: equatorial_radius_m,
            flattening,
        })
    }
    pub fn equatorial_radius_m(self) -> f64 {
        self.a
    }
    pub fn flattening(self) -> f64 {
        self.flattening
    }
    /// Surface radii at geodetic latitude; height and longitude are ignored.
    /// These are local normal-section curvatures, not global ray-path lengths.
    pub fn curvature_radii(self, point: Geodetic) -> Result<CurvatureRadii> {
        let (sin_lat, cos_lat) = point.lat.to_radians().sin_cos();
        let axis_ratio = 1.0 - self.flattening;
        // Avoid 1 - e² sin²(lat) cancellation on very flattened ellipsoids.
        let q = cos_lat.hypot(axis_ratio * sin_lat);
        let prime_vertical_m = self.a / q;
        let meridional_m = prime_vertical_m * (axis_ratio / q).powi(2);
        if !prime_vertical_m.is_finite() || !meridional_m.is_finite() || meridional_m <= 0.0 {
            return Err(Error::InvalidInput(
                "curvature exceeds numerical range".into(),
            ));
        }
        Ok(CurvatureRadii {
            meridional_m,
            prime_vertical_m,
        })
    }
    pub fn meridional_radius_m(self, point: Geodetic) -> Result<f64> {
        Ok(self.curvature_radii(point)?.meridional_m)
    }
    pub fn prime_vertical_radius_m(self, point: Geodetic) -> Result<f64> {
        Ok(self.curvature_radii(point)?.prime_vertical_m)
    }
    /// Euler normal-section radius at an azimuth clockwise from north.
    pub fn radius_along_azimuth_m(self, point: Geodetic, azimuth_deg: f64) -> Result<f64> {
        self.curvature_radii(point)?.along_azimuth_m(azimuth_deg)
    }
    pub fn to_cartesian(self, point: Geodetic) -> Ecef {
        let e2 = self.flattening * (2.0 - self.flattening);
        let (sl, cl) = point.lat.to_radians().sin_cos();
        let (so, co) = point.lon.to_radians().sin_cos();
        let n = self.a / (1.0 - e2 * sl * sl).sqrt();
        Ecef {
            x: (n + point.height) * cl * co,
            y: (n + point.height) * cl * so,
            z: (n * (1.0 - e2) + point.height) * sl,
        }
    }
    pub fn to_geodetic(self, point: Ecef) -> Result<Geodetic> {
        // Scale the general ellipsoid inversion; no WGS84 constants are used.
        finite(point.x, "x")?;
        finite(point.y, "y")?;
        finite(point.z, "z")?;
        let p = point.x.hypot(point.y);
        if p.hypot(point.z) < self.a / 2.0 {
            return Err(Error::InvalidInput(
                "coordinate too near ellipsoid centre".into(),
            ));
        }
        if p < self.a * 1e-15 {
            return Geodetic::new(
                90.0_f64.copysign(point.z),
                0.0,
                point.z.abs() - self.a * (1.0 - self.flattening),
            );
        }
        let e2 = self.flattening * (2.0 - self.flattening);
        let mut lat = point.z.atan2(p * (1.0 - e2));
        for _ in 0..30 {
            let n = self.a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
            let next = (point.z + e2 * n * lat.sin()).atan2(p);
            if (next - lat).abs() < 1e-14 {
                lat = next;
                break;
            }
            lat = next;
        }
        let n = self.a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
        let h = p * lat.cos() + point.z * lat.sin() - n * (1.0 - e2 * lat.sin().powi(2));
        Geodetic::new(lat.to_degrees(), point.y.atan2(point.x).to_degrees(), h)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvatureRadii {
    meridional_m: f64,
    prime_vertical_m: f64,
}
impl CurvatureRadii {
    /// Radius of the north/south meridian section, in metres.
    pub fn meridional_m(self) -> f64 {
        self.meridional_m
    }
    /// Radius of the east/west prime-vertical section, in metres.
    pub fn prime_vertical_m(self) -> f64 {
        self.prime_vertical_m
    }
    /// Euler's formula: 1/R = cos²(azimuth)/M + sin²(azimuth)/N.
    pub fn along_azimuth_m(self, azimuth_deg: f64) -> Result<f64> {
        finite(azimuth_deg, "curvature azimuth")?;
        let (sin_az, cos_az) = azimuth_deg.to_radians().sin_cos();
        // Scaled harmonic mean avoids reciprocal under/overflow.
        let ratio = self.meridional_m / self.prime_vertical_m;
        Ok(self.meridional_m / (cos_az * cos_az + ratio * sin_az * sin_az))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodetic {
    lat: f64,
    lon: f64,
    height: f64,
}

impl Geodetic {
    pub fn new(latitude_deg: f64, longitude_deg: f64, ellipsoidal_height_m: f64) -> Result<Self> {
        finite(latitude_deg, "latitude")?;
        finite(longitude_deg, "longitude")?;
        finite(ellipsoidal_height_m, "height")?;
        if !(-90.0..=90.0).contains(&latitude_deg) {
            return Err(Error::InvalidInput("latitude must be in [-90, 90]".into()));
        }
        Ok(Self {
            lat: latitude_deg,
            lon: if (-180.0..180.0).contains(&longitude_deg) {
                longitude_deg
            } else {
                (longitude_deg + 180.0).rem_euclid(360.0) - 180.0
            },
            height: ellipsoidal_height_m,
        })
    }
    /// Plain constructor for callers that have already validated their inputs.
    /// Validation is still enforced; invalid latitude or nonfinite input panics.
    /// Use `new` for data read from files, users or sensing feeds.
    ///
    /// # Panics
    /// Panics if the inputs do not satisfy `Geodetic::new`'s contract.
    pub fn from_validated(
        latitude_deg: f64,
        longitude_deg: f64,
        ellipsoidal_height_m: f64,
    ) -> Self {
        Self::new(latitude_deg, longitude_deg, ellipsoidal_height_m)
            .expect("prevalidated geodetic coordinate")
    }
    /// WGS84 surface curvatures for this validated latitude; height is ignored.
    pub fn curvature_radii(self) -> CurvatureRadii {
        Ellipsoid::WGS84
            .curvature_radii(self)
            .expect("finite WGS84 curvature")
    }
    pub fn latitude_deg(self) -> f64 {
        self.lat
    }
    pub fn longitude_deg(self) -> f64 {
        self.lon
    }
    pub fn height_m(self) -> f64 {
        self.height
    }
    pub fn with_height(self, height_m: f64) -> Result<Self> {
        Self::new(self.lat, self.lon, height_m)
    }

    pub fn to_ecef(self) -> Ecef {
        let (sl, cl) = self.lat.to_radians().sin_cos();
        let (so, co) = self.lon.to_radians().sin_cos();
        let n = WGS84_A / (1.0 - E2 * sl * sl).sqrt();
        Ecef {
            x: (n + self.height) * cl * co,
            y: (n + self.height) * cl * so,
            z: (n * (1.0 - E2) + self.height) * sl,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ecef {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Representation conversion only; finite-value validation remains the consumer's responsibility.
impl From<[f64; 3]> for Ecef {
    fn from([x, y, z]: [f64; 3]) -> Self {
        Self { x, y, z }
    }
}
impl From<Ecef> for [f64; 3] {
    fn from(point: Ecef) -> Self {
        [point.x, point.y, point.z]
    }
}
impl From<Vec3> for Ecef {
    fn from(vector: Vec3) -> Self {
        Self {
            x: vector.x,
            y: vector.y,
            z: vector.z,
        }
    }
}
impl From<Ecef> for Vec3 {
    fn from(point: Ecef) -> Self {
        Self {
            x: point.x,
            y: point.y,
            z: point.z,
        }
    }
}
/// Difference of positions is a displacement vector.
impl Sub for Ecef {
    type Output = Vec3;
    fn sub(self, other: Self) -> Vec3 {
        Vec3::from(self) - Vec3::from(other)
    }
}
impl Add<Vec3> for Ecef {
    type Output = Self;
    fn add(self, displacement: Vec3) -> Self {
        (Vec3::from(self) + displacement).into()
    }
}
/// Cartesian component sum, useful for weighted/chord averages. This is not a
/// surface-geodesic midpoint; use `geodesics::GeodesicPath` for that operation.
impl Add for Ecef {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        (Vec3::from(self) + Vec3::from(other)).into()
    }
}
impl Sub<Vec3> for Ecef {
    type Output = Self;
    fn sub(self, displacement: Vec3) -> Self {
        (Vec3::from(self) - displacement).into()
    }
}
impl Mul<f64> for Ecef {
    type Output = Self;
    fn mul(self, scale: f64) -> Self {
        (Vec3::from(self) * scale).into()
    }
}
impl Div<f64> for Ecef {
    type Output = Self;
    fn div(self, scale: f64) -> Self {
        Self {
            x: self.x / scale,
            y: self.y / scale,
            z: self.z / scale,
        }
    }
}

impl Ecef {
    pub fn to_geodetic(self) -> Result<Geodetic> {
        finite(self.x, "ECEF x")?;
        finite(self.y, "ECEF y")?;
        finite(self.z, "ECEF z")?;
        let p = self.x.hypot(self.y);
        // Geodetic latitude is ambiguous deep inside the ellipsoid. This API is
        // intended for surface and above-surface coordinates, not Earth's interior.
        if p.hypot(self.z) < WGS84_A / 2.0 {
            return Err(Error::InvalidInput(
                "ECEF point is too near Earth's centre".into(),
            ));
        }
        if p < 1e-8 {
            return Geodetic::new(
                90.0_f64.copysign(self.z),
                0.0,
                self.z.abs() - WGS84_A * (1.0 - WGS84_F),
            );
        }
        let mut lat = self.z.atan2(p * (1.0 - E2));
        for _ in 0..20 {
            let n = WGS84_A / (1.0 - E2 * lat.sin().powi(2)).sqrt();
            let next = (self.z + E2 * n * lat.sin()).atan2(p);
            if (next - lat).abs() < 1e-14 {
                lat = next;
                break;
            }
            lat = next;
        }
        let n = WGS84_A / (1.0 - E2 * lat.sin().powi(2)).sqrt();
        let height = p * lat.cos() + self.z * lat.sin() - n * (1.0 - E2 * lat.sin().powi(2));
        Geodetic::new(lat.to_degrees(), self.y.atan2(self.x).to_degrees(), height)
    }
    pub fn distance(self, other: Self) -> f64 {
        (self.x - other.x)
            .hypot(self.y - other.y)
            .hypot(self.z - other.z)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Enu {
    pub east: f64,
    pub north: f64,
    pub up: f64,
}

/// Tangent frame at a WGS84 origin. ENU is a local Cartesian frame, not a map projection.
pub struct LocalFrame {
    origin: Geodetic,
}
impl LocalFrame {
    pub fn new(origin: Geodetic) -> Self {
        Self { origin }
    }
    pub fn to_enu(&self, point: Ecef) -> Enu {
        let o = self.origin.to_ecef();
        let (sl, cl) = self.origin.lat.to_radians().sin_cos();
        let (so, co) = self.origin.lon.to_radians().sin_cos();
        let (x, y, z) = (point.x - o.x, point.y - o.y, point.z - o.z);
        Enu {
            east: -so * x + co * y,
            north: -sl * co * x - sl * so * y + cl * z,
            up: cl * co * x + cl * so * y + sl * z,
        }
    }
    pub fn to_ecef(&self, point: Enu) -> Ecef {
        let o = self.origin.to_ecef();
        let (sl, cl) = self.origin.lat.to_radians().sin_cos();
        let (so, co) = self.origin.lon.to_radians().sin_cos();
        Ecef {
            x: o.x - so * point.east - sl * co * point.north + cl * co * point.up,
            y: o.y + co * point.east - sl * so * point.north + cl * so * point.up,
            z: o.z + cl * point.north + sl * point.up,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hemisphere {
    North,
    South,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Utm {
    zone: u8,
    hemisphere: Hemisphere,
}

impl Utm {
    pub fn new(zone: u8, hemisphere: Hemisphere) -> Result<Self> {
        if !(1..=60).contains(&zone) {
            return Err(Error::InvalidInput("UTM zone must be 1..=60".into()));
        }
        Ok(Self { zone, hemisphere })
    }
    pub fn zone(self) -> u8 {
        self.zone
    }
    pub fn hemisphere(self) -> Hemisphere {
        self.hemisphere
    }
    fn central_meridian(self) -> f64 {
        (f64::from(self.zone) * 6.0 - 183.0).to_radians()
    }
    pub fn project(self, point: Geodetic) -> Result<(f64, f64)> {
        if !(-80.0..=84.0).contains(&point.lat)
            || (point.lat < 0.0 && self.hemisphere == Hemisphere::North)
            || (point.lat > 0.0 && self.hemisphere == Hemisphere::South)
        {
            return Err(Error::InvalidInput(
                "latitude outside selected UTM hemisphere or range".into(),
            ));
        }
        let phi = point.lat.to_radians();
        let delta = (point.lon.to_radians() - self.central_meridian() + std::f64::consts::PI)
            .rem_euclid(2.0 * std::f64::consts::PI)
            - std::f64::consts::PI;
        if delta.abs() > 6.0_f64.to_radians() {
            return Err(Error::InvalidInput(
                "point is more than 6 degrees from UTM central meridian".into(),
            ));
        }
        let ep2 = E2 / (1.0 - E2);
        let n = WGS84_A / (1.0 - E2 * phi.sin().powi(2)).sqrt();
        let t = phi.tan().powi(2);
        let c = ep2 * phi.cos().powi(2);
        let a = phi.cos() * delta;
        let m = WGS84_A
            * ((1.0 - E2 / 4.0 - 3.0 * E2.powi(2) / 64.0 - 5.0 * E2.powi(3) / 256.0) * phi
                - (3.0 * E2 / 8.0 + 3.0 * E2.powi(2) / 32.0 + 45.0 * E2.powi(3) / 1024.0)
                    * (2.0 * phi).sin()
                + (15.0 * E2.powi(2) / 256.0 + 45.0 * E2.powi(3) / 1024.0) * (4.0 * phi).sin()
                - 35.0 * E2.powi(3) / 3072.0 * (6.0 * phi).sin());
        let easting = 500_000.0
            + 0.9996
                * n
                * (a + (1.0 - t + c) * a.powi(3) / 6.0
                    + (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * ep2) * a.powi(5) / 120.0);
        let northing = 0.9996
            * (m + n
                * phi.tan()
                * (a * a / 2.0
                    + (5.0 - t + 9.0 * c + 4.0 * c * c) * a.powi(4) / 24.0
                    + (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * ep2) * a.powi(6) / 720.0))
            + if self.hemisphere == Hemisphere::South {
                10_000_000.0
            } else {
                0.0
            };
        Ok((easting, northing))
    }
    pub fn unproject(self, easting: f64, northing: f64, height_m: f64) -> Result<Geodetic> {
        finite(easting, "easting")?;
        finite(northing, "northing")?;
        if !(100_000.0..=900_000.0).contains(&easting) || !(0.0..=10_000_000.0).contains(&northing)
        {
            return Err(Error::InvalidInput(
                "UTM coordinates outside supported range".into(),
            ));
        }
        let x = easting - 500_000.0;
        let y = northing
            - if self.hemisphere == Hemisphere::South {
                10_000_000.0
            } else {
                0.0
            };
        let mu = y
            / (0.9996
                * WGS84_A
                * (1.0 - E2 / 4.0 - 3.0 * E2.powi(2) / 64.0 - 5.0 * E2.powi(3) / 256.0));
        let e1 = (1.0 - (1.0 - E2).sqrt()) / (1.0 + (1.0 - E2).sqrt());
        let fp = mu
            + (3.0 * e1 / 2.0 - 27.0 * e1.powi(3) / 32.0) * (2.0 * mu).sin()
            + (21.0 * e1.powi(2) / 16.0 - 55.0 * e1.powi(4) / 32.0) * (4.0 * mu).sin()
            + 151.0 * e1.powi(3) / 96.0 * (6.0 * mu).sin()
            + 1097.0 * e1.powi(4) / 512.0 * (8.0 * mu).sin();
        let ep2 = E2 / (1.0 - E2);
        let c = ep2 * fp.cos().powi(2);
        let t = fp.tan().powi(2);
        let n = WGS84_A / (1.0 - E2 * fp.sin().powi(2)).sqrt();
        let r = WGS84_A * (1.0 - E2) / (1.0 - E2 * fp.sin().powi(2)).powf(1.5);
        let d = x / (n * 0.9996);
        let lat = fp
            - n * fp.tan() / r
                * (d * d / 2.0
                    - (5.0 + 3.0 * t + 10.0 * c - 4.0 * c * c - 9.0 * ep2) * d.powi(4) / 24.0
                    + (61.0 + 90.0 * t + 298.0 * c + 45.0 * t * t - 252.0 * ep2 - 3.0 * c * c)
                        * d.powi(6)
                        / 720.0);
        let lon = self.central_meridian()
            + (d - (1.0 + 2.0 * t + c) * d.powi(3) / 6.0
                + (5.0 - 2.0 * c + 28.0 * t - 3.0 * c * c + 8.0 * ep2 + 24.0 * t * t) * d.powi(5)
                    / 120.0)
                / fp.cos();
        let point = Geodetic::new(lat.to_degrees(), lon.to_degrees(), height_m)?;
        self.project(point)?;
        Ok(point)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Crs {
    Wgs84,
    /// Longitude/latitude on the ellipsoid of an explicitly associated body.
    /// Angular mapping only: the body and ellipsoid live in the atlas/frame metadata.
    Planetographic,
    WebMercator,
    Utm(Utm),
}
impl Crs {
    /// WGS84 raster coordinates use (longitude, latitude), unlike Geodetic::new.
    pub fn project(self, point: Geodetic) -> Result<(f64, f64)> {
        match self {
            Self::Wgs84 | Self::Planetographic => Ok((point.lon, point.lat)),
            Self::Utm(utm) => utm.project(point),
            Self::WebMercator => {
                if point.lat.abs() > 85.051_128_779_806_6 {
                    return Err(Error::OutsideCoverage);
                }
                Ok((
                    WGS84_A * point.lon.to_radians(),
                    WGS84_A
                        * (std::f64::consts::FRAC_PI_4 + point.lat.to_radians() / 2.0)
                            .tan()
                            .ln(),
                ))
            }
        }
    }
    pub fn unproject(self, x: f64, y: f64, height_m: f64) -> Result<Geodetic> {
        finite(x, "x")?;
        finite(y, "y")?;
        match self {
            Self::Wgs84 | Self::Planetographic => Geodetic::new(y, x, height_m),
            Self::Utm(utm) => utm.unproject(x, y, height_m),
            Self::WebMercator => {
                if x.abs() > std::f64::consts::PI * WGS84_A
                    || y.abs() > std::f64::consts::PI * WGS84_A
                {
                    return Err(Error::OutsideCoverage);
                }
                Geodetic::new(
                    (2.0 * (y / WGS84_A).exp().atan() - std::f64::consts::FRAC_PI_2).to_degrees(),
                    (x / WGS84_A).to_degrees(),
                    height_m,
                )
            }
        }
    }
}
