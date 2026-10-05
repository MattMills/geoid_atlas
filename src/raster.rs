//! In-memory georeferenced DEMs and DRG imagery, with ESRI ASCII DEM ingestion.
//! Affine transforms locate pixel CENTRES. No implicit reprojection or nodata filling.
use crate::coordinates::{Crs, Geodetic};
use crate::{Error, Result, finite};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Affine {
    origin_x: f64,
    origin_y: f64,
    col_x: f64,
    col_y: f64,
    row_x: f64,
    row_y: f64,
}
impl Affine {
    pub fn new(
        origin_x: f64,
        origin_y: f64,
        col_x: f64,
        col_y: f64,
        row_x: f64,
        row_y: f64,
    ) -> Result<Self> {
        for v in [origin_x, origin_y, col_x, col_y, row_x, row_y] {
            finite(v, "affine coefficient")?;
        }
        let determinant = col_x * row_y - row_x * col_y;
        if !determinant.is_finite() || determinant == 0.0 {
            return Err(Error::InvalidInput(
                "affine transform must be invertible".into(),
            ));
        }
        Ok(Self {
            origin_x,
            origin_y,
            col_x,
            col_y,
            row_x,
            row_y,
        })
    }
    pub fn world(self, col: f64, row: f64) -> (f64, f64) {
        (
            self.origin_x + col * self.col_x + row * self.row_x,
            self.origin_y + col * self.col_y + row * self.row_y,
        )
    }
    pub fn pixel(self, x: f64, y: f64) -> Result<(f64, f64)> {
        finite(x, "x")?;
        finite(y, "y")?;
        let dx = x - self.origin_x;
        let dy = y - self.origin_y;
        let d = self.col_x * self.row_y - self.row_x * self.col_y;
        let col = (dx * self.row_y - dy * self.row_x) / d;
        let row = (dy * self.col_x - dx * self.col_y) / d;
        finite(col, "column")?;
        finite(row, "row")?;
        Ok((col, row))
    }
}

#[derive(Debug, Clone)]
pub struct Raster<T> {
    width: usize,
    height: usize,
    transform: Affine,
    crs: Crs,
    cells: Vec<Option<T>>,
}
impl<T> Raster<T> {
    pub fn new(
        width: usize,
        height: usize,
        transform: Affine,
        crs: Crs,
        cells: Vec<Option<T>>,
    ) -> Result<Self> {
        if width == 0 || height == 0 || width.checked_mul(height) != Some(cells.len()) {
            return Err(Error::InvalidInput(
                "raster dimensions must match nonempty data".into(),
            ));
        }
        Ok(Self {
            width,
            height,
            transform,
            crs,
            cells,
        })
    }
    pub fn width(&self) -> usize {
        self.width
    }
    pub fn height(&self) -> usize {
        self.height
    }
    pub fn transform(&self) -> Affine {
        self.transform
    }
    pub fn crs(&self) -> Crs {
        self.crs
    }
    pub fn cell(&self, col: usize, row: usize) -> Result<&T> {
        if col >= self.width || row >= self.height {
            return Err(Error::OutsideCoverage);
        }
        self.cells[row * self.width + col]
            .as_ref()
            .ok_or(Error::NoData)
    }
    pub fn nearest_xy(&self, x: f64, y: f64) -> Result<&T> {
        let (col, row) = self.transform.pixel(x, y)?;
        if col < -0.5
            || row < -0.5
            || col >= self.width as f64 - 0.5
            || row >= self.height as f64 - 0.5
        {
            return Err(Error::OutsideCoverage);
        }
        self.cell(col.round().max(0.0) as usize, row.round().max(0.0) as usize)
    }
    pub fn nearest(&self, point: Geodetic) -> Result<&T> {
        let (x, y) = self.crs.project(point)?;
        self.nearest_xy(x, y)
    }
}
impl Raster<f64> {
    /// Bilinear interpolation within the convex hull of sample centres.
    /// Missing cells with nonzero weight cause NoData; no extrapolation occurs.
    pub fn bilinear_xy(&self, x: f64, y: f64) -> Result<f64> {
        let (col, row) = self.transform.pixel(x, y)?;
        if col < 0.0 || row < 0.0 || col > (self.width - 1) as f64 || row > (self.height - 1) as f64
        {
            return Err(Error::OutsideCoverage);
        }
        let c0 = col.floor() as usize;
        let r0 = row.floor() as usize;
        let c1 = (c0 + 1).min(self.width - 1);
        let r1 = (r0 + 1).min(self.height - 1);
        let dx = col - c0 as f64;
        let dy = row - r0 as f64;
        let mut value = 0.0;
        for (c, r, w) in [
            (c0, r0, (1.0 - dx) * (1.0 - dy)),
            (c1, r0, dx * (1.0 - dy)),
            (c0, r1, (1.0 - dx) * dy),
            (c1, r1, dx * dy),
        ] {
            if w > 0.0 {
                let v = *self.cell(c, r)?;
                finite(v, "raster value")?;
                value += v * w;
            }
        }
        finite(value, "interpolated value")?;
        Ok(value)
    }
    pub fn bilinear(&self, point: Geodetic) -> Result<f64> {
        let (x, y) = self.crs.project(point)?;
        self.bilinear_xy(x, y)
    }
}

/// Vertical references must be selected explicitly when ingesting elevations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerticalDatum {
    Ellipsoid,
    Orthometric,
}
pub struct Dem {
    pub raster: Raster<f64>,
    pub vertical_datum: VerticalDatum,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rgb {
    pub red: u8,
    pub green: u8,
    pub blue: u8,
}
/// Already georeferenced digital raster graphic; sampling does not infer semantics.
pub type Drg = Raster<Rgb>;

impl Dem {
    /// ESRI ASCII grid, north-to-south rows; external CRS/datum are mandatory.
    /// Supports xllcorner/yllcorner or xllcenter/yllcenter, optional NODATA_value.
    pub fn from_esri_ascii(text: &str, crs: Crs, vertical_datum: VerticalDatum) -> Result<Self> {
        let mut header: HashMap<String, String> = HashMap::new();
        let mut data = Vec::new();
        let mut reading_data = false;
        for line in text.lines().filter(|s| !s.trim().is_empty()) {
            let parts: Vec<_> = line.split_whitespace().collect();
            let key = parts[0].to_ascii_lowercase();
            let known = matches!(
                key.as_str(),
                "ncols"
                    | "nrows"
                    | "xllcorner"
                    | "yllcorner"
                    | "xllcenter"
                    | "yllcenter"
                    | "cellsize"
                    | "nodata_value"
            );
            if known && !reading_data {
                if parts.len() != 2 || header.insert(key, parts[1].into()).is_some() {
                    return Err(Error::Parse("invalid or duplicate grid header".into()));
                }
            } else {
                reading_data = true;
                for p in parts {
                    data.push(
                        p.parse::<f64>()
                            .map_err(|_| Error::Parse(format!("invalid grid value: {p}")))?,
                    );
                }
            }
        }
        let number = |key: &str| -> Result<f64> {
            let v = header
                .get(key)
                .ok_or_else(|| Error::Parse(format!("missing {key}")))?
                .parse::<f64>()
                .map_err(|_| Error::Parse(format!("invalid {key}")))?;
            finite(v, key)?;
            Ok(v)
        };
        let dimension = |key: &str| -> Result<usize> {
            header
                .get(key)
                .ok_or_else(|| Error::Parse(format!("missing {key}")))?
                .parse::<usize>()
                .map_err(|_| Error::Parse(format!("invalid {key}")))
        };
        let width = dimension("ncols")?;
        let height = dimension("nrows")?;
        let size = number("cellsize")?;
        if size <= 0.0 {
            return Err(Error::Parse("cellsize must be positive".into()));
        }
        let lower = |corner: &str, center: &str| -> Result<f64> {
            match (header.contains_key(corner), header.contains_key(center)) {
                (true, false) => Ok(number(corner)? + size / 2.0),
                (false, true) => number(center),
                _ => Err(Error::Parse(format!(
                    "supply exactly one of {corner} or {center}"
                ))),
            }
        };
        let x = lower("xllcorner", "xllcenter")?;
        let y = lower("yllcorner", "yllcenter")?;
        let nodata = header
            .get("nodata_value")
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| Error::Parse("invalid nodata value".into()))
            })
            .transpose()?;
        let cells = data
            .into_iter()
            .map(|v| {
                if nodata.is_some_and(|n| n == v || (n.is_nan() && v.is_nan())) {
                    Ok(None)
                } else {
                    finite(v, "DEM elevation")?;
                    Ok(Some(v))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        let transform = Affine::new(
            x,
            y + height.saturating_sub(1) as f64 * size,
            size,
            0.0,
            0.0,
            -size,
        )?;
        Ok(Self {
            raster: Raster::new(width, height, transform, crs, cells)?,
            vertical_datum,
        })
    }
}
