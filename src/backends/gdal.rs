//! GDAL >= 3.8 local raster ingestion and strict GeoTIFF warping.
//! Metadata retains the native CRS; it is never relabeled as the core WGS84 CRS.
use super::{Handle, backend_error, copy_string, cstring};
use crate::raster::Affine;
use crate::{Error, Result, finite};
use std::ffi::{c_char, c_int, c_void};
use std::path::{Path, PathBuf};
use std::ptr;
use std::sync::Once;

unsafe extern "C" {
    fn GDALAllRegister();
    fn GDALVersionInfo(request: *const c_char) -> *const c_char;
    fn GDALOpen(path: *const c_char, access: c_int) -> *mut c_void;
    fn GDALClose(dataset: *mut c_void) -> c_int;
    fn GDALGetRasterXSize(dataset: *mut c_void) -> c_int;
    fn GDALGetRasterYSize(dataset: *mut c_void) -> c_int;
    fn GDALGetRasterCount(dataset: *mut c_void) -> c_int;
    fn GDALGetProjectionRef(dataset: *mut c_void) -> *const c_char;
    fn GDALGetGeoTransform(dataset: *mut c_void, transform: *mut f64) -> c_int;
    fn GDALGetRasterBand(dataset: *mut c_void, band: c_int) -> *mut c_void;
    fn GDALGetRasterNoDataValue(band: *mut c_void, success: *mut c_int) -> f64;
    fn GDALGetRasterScale(band: *mut c_void, success: *mut c_int) -> f64;
    fn GDALGetRasterOffset(band: *mut c_void, success: *mut c_int) -> f64;
    fn GDALGetRasterUnitType(band: *mut c_void) -> *const c_char;
    fn GDALGetMaskFlags(band: *mut c_void) -> c_int;
    fn GDALGetMaskBand(band: *mut c_void) -> *mut c_void;
    fn GDALRasterIO(
        band: *mut c_void,
        mode: c_int,
        x: c_int,
        y: c_int,
        width: c_int,
        height: c_int,
        buffer: *mut c_void,
        buffer_width: c_int,
        buffer_height: c_int,
        data_type: c_int,
        pixel_space: c_int,
        line_space: c_int,
    ) -> c_int;
    fn CPLGetLastErrorMsg() -> *const c_char;
    fn CPLErrorReset();
    fn OSRGetPROJEnableNetwork() -> c_int;
    fn OSRGetPROJVersion(major: *mut c_int, minor: *mut c_int, patch: *mut c_int);
    fn GDALWarpAppOptionsNew(argv: *mut *mut c_char, binary_options: *mut c_void) -> *mut c_void;
    fn GDALWarpAppOptionsFree(options: *mut c_void);
    fn GDALWarp(
        destination: *const c_char,
        output: *mut c_void,
        count: c_int,
        sources: *mut *mut c_void,
        options: *const c_void,
        usage_error: *mut c_int,
    ) -> *mut c_void;
}
unsafe extern "C" fn close(dataset: *mut c_void) {
    // SAFETY: Handle owns this live dataset and closes it once.
    unsafe {
        GDALClose(dataset);
    }
}
fn failure() -> Error {
    // SAFETY: GDAL thread-local error string is copied immediately.
    unsafe { backend_error("GDAL", copy_string(CPLGetLastErrorMsg())) }
}
fn initialize() -> Result<()> {
    let request = cstring("VERSION_NUM")?;
    // SAFETY: Static version string is copied; Once serializes driver registration.
    unsafe {
        let version = copy_string(GDALVersionInfo(request.as_ptr()))
            .parse::<u32>()
            .unwrap_or(0);
        if version < 3_080_000 {
            return Err(backend_error("GDAL", "GDAL >= 3.8 required"));
        }
        static REGISTER: Once = Once::new();
        REGISTER.call_once(|| GDALAllRegister());
    }
    Ok(())
}
pub fn version() -> String {
    // SAFETY: Request is static NUL-terminated, result copied.
    unsafe { copy_string(GDALVersionInfo(c"RELEASE_NAME".as_ptr())) }
}
#[derive(Debug, Clone)]
pub struct RasterMetadata {
    pub width: usize,
    pub height: usize,
    pub bands: usize,
    pub crs_wkt: String,
    /// GDAL corner-based affine coefficients. A window exposes pixel CENTRES.
    pub corner_transform: [f64; 6],
}
pub struct GdalRaster {
    dataset: Handle,
    metadata: RasterMetadata,
}
#[derive(Debug, Clone)]
pub struct RasterWindow {
    pub width: usize,
    pub height: usize,
    pub crs_wkt: String,
    /// Native band unit label, possibly empty/unspecified. Never inferred from the CRS.
    pub unit: String,
    pub transform: Affine,
    /// Row-major physical values (native value * scale + offset). Mask/nodata/nonfinite -> None.
    pub cells: Vec<Option<f64>>,
}
impl RasterWindow {
    /// Connect a window to the core raster/atlas only after native CRS equivalence
    /// validation. This never reprojects cells. Warp first if references differ.
    #[cfg(feature = "proj-backend")]
    pub fn into_core_raster(
        self,
        crs: crate::coordinates::Crs,
    ) -> Result<crate::raster::Raster<f64>> {
        use crate::coordinates::{Crs, Hemisphere};
        let target = match crs {
            Crs::Wgs84 => "EPSG:4326".to_string(),
            Crs::WebMercator => "EPSG:3857".to_string(),
            Crs::Utm(utm) => format!(
                "EPSG:{}",
                match utm.hemisphere() {
                    Hemisphere::North => 32600,
                    Hemisphere::South => 32700,
                } + u32::from(utm.zone())
            ),
            Crs::Planetographic => {
                return Err(Error::InvalidInput(
                    "planetary CRS needs an explicit body adapter".into(),
                ));
            }
        };
        if !super::proj::equivalent_crs(&self.crs_wkt, &target)? {
            return Err(Error::InvalidInput(
                "raster CRS differs; warp before attaching to core atlas".into(),
            ));
        }
        crate::raster::Raster::new(self.width, self.height, self.transform, crs, self.cells)
    }
    pub fn nearest(&self, x: f64, y: f64) -> Result<f64> {
        let (col, row) = self.transform.pixel(x, y)?;
        let col = col.round();
        let row = row.round();
        if col < 0.0 || row < 0.0 || col >= self.width as f64 || row >= self.height as f64 {
            return Err(Error::OutsideCoverage);
        }
        self.cells[row as usize * self.width + col as usize].ok_or(Error::NoData)
    }
}
impl GdalRaster {
    /// Opens a local existing file read-only. Retains driver/CRS behavior from GDAL.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        initialize()?;
        let path = path.as_ref();
        if !path.is_file() {
            return Err(backend_error(
                "GDAL",
                "input must be an existing local file",
            ));
        }
        let path = cstring(
            path.to_str()
                .ok_or_else(|| Error::InvalidInput("raster path must be UTF-8".into()))?,
        )?;
        // SAFETY: Dataset is owned by Handle; native metadata strings are copied.
        unsafe {
            CPLErrorReset();
            let dataset = Handle::new(
                GDALOpen(path.as_ptr(), 0),
                close,
                "GDAL",
                "cannot open dataset",
            )
            .map_err(|_| failure())?;
            let width = GDALGetRasterXSize(dataset.ptr());
            let height = GDALGetRasterYSize(dataset.ptr());
            let bands = GDALGetRasterCount(dataset.ptr());
            if width <= 0 || height <= 0 || bands <= 0 {
                return Err(backend_error("GDAL", "empty/non-raster dataset"));
            }
            let mut corner_transform = [0.0; 6];
            if GDALGetGeoTransform(dataset.ptr(), corner_transform.as_mut_ptr()) != 0 {
                return Err(failure());
            }
            for value in corner_transform {
                finite(value, "raster geotransform")?;
            }
            let crs_wkt = copy_string(GDALGetProjectionRef(dataset.ptr()));
            if crs_wkt.is_empty() {
                return Err(backend_error("GDAL", "raster has no CRS"));
            }
            let metadata = RasterMetadata {
                width: width as usize,
                height: height as usize,
                bands: bands as usize,
                corner_transform,
                crs_wkt,
            };
            Ok(Self { dataset, metadata })
        }
    }
    pub fn metadata(&self) -> &RasterMetadata {
        &self.metadata
    }
    /// Read one 1-based band, bounded to the source window. No implicit resampling.
    pub fn read_window(
        &mut self,
        band_index: usize,
        col: usize,
        row: usize,
        width: usize,
        height: usize,
    ) -> Result<RasterWindow> {
        let m = &self.metadata;
        if band_index == 0
            || band_index > m.bands
            || width == 0
            || height == 0
            || col.checked_add(width).is_none_or(|end| end > m.width)
            || row.checked_add(height).is_none_or(|end| end > m.height)
        {
            return Err(Error::OutsideCoverage);
        }
        let count = width
            .checked_mul(height)
            .ok_or_else(|| Error::InvalidInput("window overflow".into()))?;
        let mut values = vec![0.0_f64; count];
        let mut mask = vec![255_u8; count];
        // SAFETY: Bounds fit GDAL's signed source dimensions; buffers have width*height
        // elements of the requested Float64/Byte types. Bands are borrowed from dataset.
        unsafe {
            let band = GDALGetRasterBand(self.dataset.ptr(), band_index as c_int);
            if band.is_null() {
                return Err(failure());
            }
            if GDALRasterIO(
                band,
                0,
                col as c_int,
                row as c_int,
                width as c_int,
                height as c_int,
                values.as_mut_ptr().cast(),
                width as c_int,
                height as c_int,
                7,
                0,
                0,
            ) != 0
            {
                return Err(failure());
            }
            if GDALGetMaskFlags(band) & 1 == 0 {
                let mask_band = GDALGetMaskBand(band);
                if mask_band.is_null()
                    || GDALRasterIO(
                        mask_band,
                        0,
                        col as c_int,
                        row as c_int,
                        width as c_int,
                        height as c_int,
                        mask.as_mut_ptr().cast(),
                        width as c_int,
                        height as c_int,
                        1,
                        0,
                        0,
                    ) != 0
                {
                    return Err(failure());
                }
            }
            let mut valid = 0;
            let nodata = GDALGetRasterNoDataValue(band, &mut valid);
            let has_nodata = valid != 0;
            let mut scale = GDALGetRasterScale(band, &mut valid);
            if valid == 0 {
                scale = 1.0;
            }
            let mut offset = GDALGetRasterOffset(band, &mut valid);
            if valid == 0 {
                offset = 0.0;
            }
            finite(scale, "raster scale")?;
            finite(offset, "raster offset")?;
            let cells = values
                .into_iter()
                .zip(mask)
                .map(|(value, mask)| {
                    if mask == 0 || !value.is_finite() || (has_nodata && value == nodata) {
                        None
                    } else {
                        let value = value * scale + offset;
                        value.is_finite().then_some(value)
                    }
                })
                .collect();
            let [x, a, b, y, d, e] = m.corner_transform;
            let transform = Affine::new(
                x + (col as f64 + 0.5) * a + (row as f64 + 0.5) * b,
                y + (col as f64 + 0.5) * d + (row as f64 + 0.5) * e,
                a,
                d,
                b,
                e,
            )?;
            Ok(RasterWindow {
                width,
                height,
                crs_wkt: m.crs_wkt.clone(),
                unit: copy_string(GDALGetRasterUnitType(band)),
                transform,
                cells,
            })
        }
    }
}
#[derive(Debug, Clone, Copy)]
pub enum Resampling {
    Nearest,
    Bilinear,
    Cubic,
    Average,
}
impl Resampling {
    fn name(self) -> &'static str {
        match self {
            Self::Nearest => "near",
            Self::Bilinear => "bilinear",
            Self::Cubic => "cubic",
            Self::Average => "average",
        }
    }
}
pub struct WarpRequest<'a> {
    pub target_crs: &'a str,
    /// Target CRS units/order: xmin, ymin, xmax, ymax. No antimeridian wrapping.
    pub bounds: [f64; 4],
    pub width: usize,
    pub height: usize,
    pub resampling: Resampling,
    pub source_decimal_year: Option<f64>,
    pub target_decimal_year: Option<f64>,
}
/// Warp into a NEW directory as raster.tif; never updates an existing destination.
/// On failure the new directory may contain partial output; it is not deleted.
/// Uses GDAL's source georeferencing, exact transformer (-et 0), offline PROJ,
/// ALLOW_BALLPARK=NO and ONLY_BEST=YES. Native PROJ >= 9.2 must be installed.
/// This is georeferenced raster warping, not raw-camera orthorectification.
pub fn warp_to_new_directory(
    source: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    request: WarpRequest<'_>,
) -> Result<PathBuf> {
    for value in request.bounds {
        finite(value, "warp bound")?;
    }
    if request.bounds[0] >= request.bounds[2]
        || request.bounds[1] >= request.bounds[3]
        || request.width == 0
        || request.height == 0
        || request.width > c_int::MAX as usize
        || request.height > c_int::MAX as usize
    {
        return Err(Error::InvalidInput(
            "ordered bounds and positive signed dimensions required".into(),
        ));
    }
    let input = GdalRaster::open(source)?;
    let mut arguments = vec![
        "-of".into(),
        "GTiff".into(),
        "-t_srs".into(),
        request.target_crs.into(),
        "-r".into(),
        request.resampling.name().into(),
        "-et".into(),
        "0".into(),
        "-to".into(),
        "ALLOW_BALLPARK=NO".into(),
        "-to".into(),
        "ONLY_BEST=YES".into(),
        "-te".into(),
    ];
    arguments.extend(request.bounds.map(|v| v.to_string()));
    arguments.extend([
        "-ts".into(),
        request.width.to_string(),
        request.height.to_string(),
    ]);
    for (name, epoch) in [
        ("-s_coord_epoch", request.source_decimal_year),
        ("-t_coord_epoch", request.target_decimal_year),
    ] {
        if let Some(epoch) = epoch {
            finite(epoch, "coordinate epoch")?;
            arguments.extend([name.into(), epoch.to_string()]);
        }
    }
    // GDAL's PROJ network setting is GLOBAL. Never mutate another application's
    // configuration; reject online settings before creating any warp transform.
    // SAFETY: Valid stack outputs; these are native configuration/version getters.
    unsafe {
        let (mut major, mut minor, mut patch) = (0, 0, 0);
        OSRGetPROJVersion(&mut major, &mut minor, &mut patch);
        if (major, minor) < (9, 2) {
            return Err(backend_error("GDAL", "warping requires PROJ >= 9.2"));
        }
        if OSRGetPROJEnableNetwork() != 0 {
            return Err(backend_error(
                "GDAL",
                "disable PROJ networking in the application before warping",
            ));
        }
    }
    let strings = arguments
        .iter()
        .map(|s| cstring(s))
        .collect::<Result<Vec<_>>>()?;
    let mut pointers: Vec<_> = strings.iter().map(|s| s.as_ptr().cast_mut()).collect();
    pointers.push(ptr::null_mut());
    let output = destination.as_ref().join("raster.tif");
    let path = cstring(
        output
            .to_str()
            .ok_or_else(|| Error::InvalidInput("output path must be UTF-8".into()))?,
    )?;
    // SAFETY: Options array is live and null-terminated; source dataset outlives warp.
    // GDAL options/output are owned by RAII handles and closed before returning.
    unsafe {
        let options = Handle::new(
            GDALWarpAppOptionsNew(pointers.as_mut_ptr(), ptr::null_mut()),
            GDALWarpAppOptionsFree,
            "GDAL",
            "invalid warp options",
        )
        .map_err(|_| failure())?;
        std::fs::create_dir(destination).map_err(|e| backend_error("GDAL", e))?;
        let mut source = input.dataset.ptr();
        let mut usage_error = 0;
        let warped = Handle::new(
            GDALWarp(
                path.as_ptr(),
                ptr::null_mut(),
                1,
                &mut source,
                options.ptr(),
                &mut usage_error,
            ),
            close,
            "GDAL",
            "warp failed",
        )
        .map_err(|_| failure())?;
        if usage_error != 0 {
            return Err(backend_error("GDAL", "invalid warp request"));
        }
        // Flush errors matter: explicitly close once and remove Handle ownership.
        let pointer = warped.ptr();
        std::mem::forget(warped);
        if GDALClose(pointer) != 0 {
            return Err(failure());
        }
    }
    Ok(output)
}
