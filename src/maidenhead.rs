//! Maidenhead grid cells, their geodetic centres and positional support.
//!
//! Standard 2-, 4-, 6- and 8-character locators are supported, plus the extended
//! 10- and 12-character convention alternating base-24 letters and base-10 digits.
//! Parsing is ASCII, case-insensitive and strict (no surrounding whitespace).
//! Maidenhead encodes longitude/latitude, not height or a survey datum: these
//! coordinates conventionally use WGS84. A centre is a representative point,
//! not a measurement accurate to the cell's centre.
//!
//! ```
//! use geoid_atlas::maidenhead::MaidenheadCell;
//! let cell: MaidenheadCell = "IO91wm".parse()?;
//! let centre = cell.centre();
//! assert_eq!(cell.locator(), "IO91WM");
//! assert!((centre.longitude_deg() + 0.125).abs() < 1e-12);
//! let uncertainty = cell.uniform_uncertainty();
//! assert!(uncertainty.east_std_m > 0.0);
//! # Ok::<(), geoid_atlas::Error>(())
//! ```

use crate::coordinates::{Ellipsoid, Geodetic};
use crate::{Error, Result};
use std::str::FromStr;

/// Geographic cell support, with bounds in degrees. Cells are west/south
/// inclusive and east/north exclusive, except the world's outer +180/+90 edges.
/// The +180 edge is equivalent to -180; retain these unwrapped bounds for display.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellBounds {
    pub west_deg: f64,
    pub south_deg: f64,
    pub east_deg: f64,
    pub north_deg: f64,
}
impl CellBounds {
    pub fn longitude_width_deg(self) -> f64 {
        self.east_deg - self.west_deg
    }
    pub fn latitude_width_deg(self) -> f64 {
        self.north_deg - self.south_deg
    }
}

/// Uncertainty **conditional on a uniform longitude/latitude position in the
/// cell**. The angular widths/bounds are exact cell support. Metre widths and
/// standard deviations are a local tangent approximation at its centre:
/// east = N cos(latitude) delta_longitude, north = M delta_latitude;
/// standard deviation = width / sqrt(12). They do not include survey, terrain,
/// datum or height error. Coarse/polar cells require a nonlinear cell model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformCellUncertainty {
    pub longitude_std_deg: f64,
    pub latitude_std_deg: f64,
    pub east_width_m: f64,
    pub north_width_m: f64,
    pub east_std_m: f64,
    pub north_std_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MaidenheadCell {
    locator: String,
    bounds: CellBounds,
}
impl MaidenheadCell {
    pub fn parse(locator: &str) -> Result<Self> {
        if !locator.is_ascii()
            || !(2..=12).contains(&locator.len())
            || !locator.len().is_multiple_of(2)
        {
            return Err(Error::Parse(
                "Maidenhead requires 2, 4, 6, 8, 10 or 12 ASCII characters".into(),
            ));
        }
        let canonical = locator.to_ascii_uppercase();
        let mut west = -180.0;
        let mut south = -90.0;
        let mut width = 360.0;
        let mut height = 180.0;
        for (level, pair) in canonical.as_bytes().as_chunks::<2>().0.iter().enumerate() {
            let (base, radix) = if level == 0 {
                (b'A', 18)
            } else if level % 2 == 1 {
                (b'0', 10)
            } else {
                (b'A', 24)
            };
            let digit = |byte: u8| -> Result<u8> {
                byte.checked_sub(base)
                    .filter(|value| *value < radix)
                    .ok_or_else(|| Error::Parse(format!("invalid Maidenhead pair {}", level + 1)))
            };
            let x = digit(pair[0])?;
            let y = digit(pair[1])?;
            width /= f64::from(radix);
            height /= f64::from(radix);
            west += f64::from(x) * width;
            south += f64::from(y) * height;
        }
        Ok(Self {
            locator: canonical,
            bounds: CellBounds {
                west_deg: west.max(-180.0),
                south_deg: south.max(-90.0),
                east_deg: (west + width).min(180.0),
                north_deg: (south + height).min(90.0),
            },
        })
    }

    pub fn locator(&self) -> &str {
        &self.locator
    }
    pub fn bounds(&self) -> CellBounds {
        self.bounds
    }
    /// Centre with height zero as a reference-surface placeholder. The locator
    /// does not tell us the site's altitude; supply a surveyed height separately.
    pub fn centre(&self) -> Geodetic {
        Geodetic::from_validated(
            (self.bounds.south_deg + self.bounds.north_deg) / 2.0,
            (self.bounds.west_deg + self.bounds.east_deg) / 2.0,
            0.0,
        )
    }
    /// Centre with an explicitly supplied **ellipsoidal** height, in metres.
    pub fn centre_at_height(&self, ellipsoidal_height_m: f64) -> Result<Geodetic> {
        self.centre().with_height(ellipsoidal_height_m)
    }
    /// WGS84 uniform-cell tangent uncertainty; see `UniformCellUncertainty`.
    pub fn uniform_uncertainty(&self) -> UniformCellUncertainty {
        self.uniform_uncertainty_on(Ellipsoid::WGS84)
            .expect("finite WGS84 cell uncertainty")
    }
    /// Same angular cell interpreted on a caller-selected ellipsoid. This changes
    /// its metre scale only; it does not perform a datum transformation.
    pub fn uniform_uncertainty_on(&self, ellipsoid: Ellipsoid) -> Result<UniformCellUncertainty> {
        let centre = self.centre();
        let radii = ellipsoid.curvature_radii(centre)?;
        let longitude_width = self.bounds.longitude_width_deg();
        let latitude_width = self.bounds.latitude_width_deg();
        let east_width_m = radii.prime_vertical_m()
            * centre.latitude_deg().to_radians().cos()
            * longitude_width.to_radians();
        let north_width_m = radii.meridional_m() * latitude_width.to_radians();
        let scale = 12.0_f64.sqrt();
        if !east_width_m.is_finite() || !north_width_m.is_finite() {
            return Err(Error::InvalidInput(
                "cell dimensions exceed numerical range".into(),
            ));
        }
        Ok(UniformCellUncertainty {
            longitude_std_deg: longitude_width / scale,
            latitude_std_deg: latitude_width / scale,
            east_width_m,
            north_width_m,
            east_std_m: east_width_m / scale,
            north_std_m: north_width_m / scale,
        })
    }
}
impl FromStr for MaidenheadCell {
    type Err = Error;
    fn from_str(locator: &str) -> Result<Self> {
        Self::parse(locator)
    }
}
