//! Strict adapters for Cartesian PROJ Helmert definitions and raster world files.
//! This is NOT a PROJ pipeline, EPSG database, WKT, or grid-shift implementation.
use crate::estimation::{Covariance6, Matrix6};
use crate::frames::{State, Vec3};
use crate::raster::Affine;
use crate::{Error, Result, finite};
use std::collections::HashMap;

type Matrix3 = [[f64; 3]; 3];
const YEAR_S: f64 = 31_557_600.0;
const ARCSECOND_RAD: f64 = std::f64::consts::PI / 648_000.0;
fn vector(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn apply(matrix: Matrix3, v: Vec3) -> Vec3 {
    let v = vector(v);
    let r = matrix.map(|row| row.iter().zip(v).map(|(a, b)| a * b).sum());
    Vec3 {
        x: r[0],
        y: r[1],
        z: r[2],
    }
}
fn inverse(m: Matrix3) -> Result<Matrix3> {
    let a = Vec3::new(m[0][0], m[0][1], m[0][2])?;
    let b = Vec3::new(m[1][0], m[1][1], m[1][2])?;
    let c = Vec3::new(m[2][0], m[2][1], m[2][2])?;
    let det = a.dot(b.cross(c));
    finite(det, "transform determinant")?;
    if det.abs() < 1e-14 {
        return Err(Error::InvalidInput("singular Cartesian transform".into()));
    }
    let columns = [
        vector(b.cross(c) * (1.0 / det)),
        vector(c.cross(a) * (1.0 / det)),
        vector(a.cross(b) * (1.0 / det)),
    ];
    let result = std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i]));
    for row in result {
        for entry in row {
            finite(entry, "inverse transform")?;
        }
    }
    Ok(result)
}

/// Instantaneous Cartesian transform y = offset + M x, with time derivative.
/// Units are metres, seconds. Velocity includes translation, scale and rotation rates.
#[derive(Debug, Clone, Copy)]
pub struct CartesianTransform {
    matrix: Matrix3,
    rate_per_s: Matrix3,
    offset: State,
}
impl CartesianTransform {
    pub fn apply(self, state: State) -> Result<State> {
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        let result = State {
            position_m: self.offset.position_m + apply(self.matrix, state.position_m),
            velocity_m_s: self.offset.velocity_m_s
                + apply(self.matrix, state.velocity_m_s)
                + apply(self.rate_per_s, state.position_m),
        };
        result.position_m.validate()?;
        result.velocity_m_s.validate()?;
        Ok(result)
    }
    pub fn unapply(self, state: State) -> Result<State> {
        state.position_m.validate()?;
        state.velocity_m_s.validate()?;
        let inv = inverse(self.matrix)?;
        let position_m = apply(inv, state.position_m - self.offset.position_m);
        let velocity_m_s = apply(
            inv,
            state.velocity_m_s - self.offset.velocity_m_s - apply(self.rate_per_s, position_m),
        );
        position_m.validate()?;
        velocity_m_s.validate()?;
        Ok(State {
            position_m,
            velocity_m_s,
        })
    }
    pub fn jacobian(self) -> Matrix6 {
        let mut result = [[0.0; 6]; 6];
        for (i, row) in result.iter_mut().take(3).enumerate() {
            row[..3].copy_from_slice(&self.matrix[i]);
        }
        for (i, row) in result.iter_mut().skip(3).enumerate() {
            row[..3].copy_from_slice(&self.rate_per_s[i]);
            row[3..].copy_from_slice(&self.matrix[i]);
        }
        result
    }
    /// Conditional covariance: transform parameter uncertainty must be added separately.
    pub fn covariance(self, covariance: Covariance6) -> Result<Covariance6> {
        covariance.propagated(self.jacobian())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RotationConvention {
    PositionVector,
    CoordinateFrame,
}
/// Linearized (small-angle) Helmert transform, matching the non-`+exact` Cartesian
/// PROJ convention. SI inside; input rotation is arcseconds, scale ppm, rates/year.
/// Supply decimal-year epochs consistent with the parameter publication. Velocity
/// rates here use a Julian year of 365.25 days; no UTC or decimal-year calendar conversion.
#[derive(Debug, Clone)]
pub struct Helmert {
    translation: Vec3,
    translation_rate: Vec3,
    rotation: Vec3,
    rotation_rate: Vec3,
    scale: f64,
    scale_rate: f64,
    reference_year: Option<f64>,
    convention: RotationConvention,
}
impl Helmert {
    /// Parse one whitespace-separated `+proj=helmert` definition. Missing numeric
    /// parameters default to zero. Unknown, duplicate and non-finite tokens fail.
    /// `+exact`, `+inv`, pipelines, grids and axis/unit remapping are unsupported.
    pub fn from_proj(text: &str) -> Result<Self> {
        let allowed = [
            "proj",
            "x",
            "y",
            "z",
            "rx",
            "ry",
            "rz",
            "s",
            "dx",
            "dy",
            "dz",
            "drx",
            "dry",
            "drz",
            "ds",
            "t_epoch",
            "convention",
        ];
        let mut tokens = HashMap::new();
        for token in text.split_whitespace() {
            let (key, value) = token
                .strip_prefix('+')
                .and_then(|t| t.split_once('='))
                .ok_or_else(|| Error::Parse("expected +key=value Helmert token".into()))?;
            if !allowed.contains(&key) || tokens.insert(key, value).is_some() {
                return Err(Error::Parse(format!(
                    "unknown or duplicate Helmert option {key}"
                )));
            }
        }
        if tokens.get("proj") != Some(&"helmert") {
            return Err(Error::Parse("expected +proj=helmert".into()));
        }
        let get = |key| -> Result<f64> {
            let value = tokens.get(key).map_or(Ok(0.0), |v| {
                v.parse::<f64>()
                    .map_err(|_| Error::Parse(format!("invalid {key}")))
            })?;
            finite(value, key)?;
            Ok(value)
        };
        let translation = Vec3::new(get("x")?, get("y")?, get("z")?)?;
        let translation_rate = Vec3::new(get("dx")?, get("dy")?, get("dz")?)?;
        let rotation = Vec3::new(get("rx")?, get("ry")?, get("rz")?)? * ARCSECOND_RAD;
        let rotation_rate = Vec3::new(get("drx")?, get("dry")?, get("drz")?)? * ARCSECOND_RAD;
        let convention = match tokens.get("convention").copied() {
            Some("position_vector") => RotationConvention::PositionVector,
            Some("coordinate_frame") => RotationConvention::CoordinateFrame,
            None if rotation == Vec3::ZERO && rotation_rate == Vec3::ZERO => {
                RotationConvention::PositionVector
            }
            _ => {
                return Err(Error::Parse(
                    "rotations require position_vector or coordinate_frame convention".into(),
                ));
            }
        };
        let result = Self {
            translation,
            translation_rate,
            rotation,
            rotation_rate,
            scale: get("s")? * 1e-6,
            scale_rate: get("ds")? * 1e-6,
            reference_year: tokens
                .contains_key("t_epoch")
                .then(|| get("t_epoch"))
                .transpose()?,
            convention,
        };
        if result.reference_year.is_none()
            && (translation_rate != Vec3::ZERO
                || rotation_rate != Vec3::ZERO
                || result.scale_rate != 0.0)
        {
            return Err(Error::Parse("rates require +t_epoch".into()));
        }
        Ok(result)
    }
    pub fn at_decimal_year(&self, year: f64) -> Result<CartesianTransform> {
        finite(year, "coordinate epoch")?;
        let dt = self
            .reference_year
            .map_or(0.0, |reference| year - reference);
        let sign = if self.convention == RotationConvention::PositionVector {
            1.0
        } else {
            -1.0
        };
        let r = (self.rotation + self.rotation_rate * dt) * sign;
        let dr = self.rotation_rate * (sign / YEAR_S);
        // Linearized Helmert formula only: large rotations require exact rotations.
        if r.norm() > 1e-3 {
            return Err(Error::InvalidInput(
                "Helmert rotation exceeds small-angle domain".into(),
            ));
        }
        let s = 1.0 + self.scale + self.scale_rate * dt;
        let ds = self.scale_rate / YEAR_S;
        finite(s, "Helmert scale")?;
        if s <= 0.0 {
            return Err(Error::InvalidInput("Helmert scale must be positive".into()));
        }
        let base = [[1.0, -r.z, r.y], [r.z, 1.0, -r.x], [-r.y, r.x, 1.0]];
        let rate = [[0.0, -dr.z, dr.y], [dr.z, 0.0, -dr.x], [-dr.y, dr.x, 0.0]];
        let matrix = base.map(|row| row.map(|x| s * x));
        let rate_per_s =
            std::array::from_fn(|i| std::array::from_fn(|j| ds * base[i][j] + s * rate[i][j]));
        let offset = State {
            position_m: self.translation + self.translation_rate * dt,
            velocity_m_s: self.translation_rate * (1.0 / YEAR_S),
        };
        offset.position_m.validate()?;
        offset.velocity_m_s.validate()?;
        for entry in matrix.into_iter().chain(rate_per_s).flatten() {
            finite(entry, "Helmert coefficient")?;
        }
        inverse(matrix)?;
        Ok(CartesianTransform {
            matrix,
            rate_per_s,
            offset,
        })
    }
}

/// Read six-line ESRI world-file coefficients A,D,B,E,C,F (pixel CENTRES).
/// A world file contains no CRS, vertical datum, units or image data: supply those separately.
pub fn read_world_file(text: &str) -> Result<Affine> {
    let values: Vec<f64> = text
        .lines()
        .map(|line| {
            line.trim()
                .parse::<f64>()
                .map_err(|_| Error::Parse("invalid world-file line".into()))
        })
        .collect::<Result<_>>()?;
    if values.len() != 6 {
        return Err(Error::Parse("world file requires six lines".into()));
    }
    Affine::new(
        values[4], values[5], values[0], values[1], values[2], values[3],
    )
}
