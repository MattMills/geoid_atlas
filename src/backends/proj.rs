//! PROJ >= 9.2 CRS transformations. Network grids and ballpark fallback are disabled.
//! Geographic x/y are longitude/latitude in degrees; other axes retain CRS units.
//! z is explicit, and t is a caller-supplied decimal coordinate year, not Epoch seconds.
//!
//! ```
//! use geoid_atlas::backends::proj::*;
//! let mut transform = ProjTransform::new("EPSG:4326", "EPSG:32613", TransformOptions::default())?;
//! let result = transform.transform(Coordinate4 {
//!     x: -105.0, y: 40.0, z: 123.0, decimal_year: 2026.0,
//! }, Direction::Forward)?;
//! assert!((result.coordinate.x - 500_000.0).abs() < 1e-6);
//! assert!(!result.operation.definition.is_empty());
//! # Ok::<(), geoid_atlas::Error>(())
//! ```
use super::{Handle, backend_error, copy_string, cstring};
use crate::{Error, Result, finite};
use std::ffi::{c_char, c_int, c_void};
use std::ptr;

#[repr(C)]
union Coordinate {
    values: [f64; 4],
}
#[repr(C)]
struct Info {
    major: c_int,
    minor: c_int,
    patch: c_int,
    release: *const c_char,
    version: *const c_char,
    searchpath: *const c_char,
    paths: *const *const c_char,
    path_count: usize,
}
unsafe extern "C" {
    fn proj_info() -> Info;
    fn proj_context_create() -> *mut c_void;
    fn proj_context_destroy(p: *mut c_void) -> *mut c_void;
    fn proj_context_set_enable_network(ctx: *mut c_void, enabled: c_int) -> c_int;
    fn proj_context_set_search_paths(ctx: *mut c_void, count: c_int, paths: *const *const c_char);
    fn proj_context_errno(ctx: *mut c_void) -> c_int;
    fn proj_context_errno_string(ctx: *mut c_void, error: c_int) -> *const c_char;
    fn proj_create(ctx: *mut c_void, definition: *const c_char) -> *mut c_void;
    fn proj_destroy(p: *mut c_void) -> *mut c_void;
    fn proj_create_crs_to_crs_from_pj(
        ctx: *mut c_void,
        source: *const c_void,
        target: *const c_void,
        area: *mut c_void,
        options: *const *const c_char,
    ) -> *mut c_void;
    fn proj_normalize_for_visualization(ctx: *mut c_void, operation: *const c_void) -> *mut c_void;
    fn proj_area_create() -> *mut c_void;
    fn proj_area_destroy(p: *mut c_void);
    fn proj_area_set_bbox(area: *mut c_void, west: f64, south: f64, east: f64, north: f64);
    fn proj_trans(p: *mut c_void, direction: c_int, coordinate: Coordinate) -> Coordinate;
    fn proj_errno(p: *const c_void) -> c_int;
    fn proj_errno_reset(p: *const c_void) -> c_int;
    fn proj_trans_get_last_used_operation(p: *mut c_void) -> *mut c_void;
    fn proj_get_name(p: *const c_void) -> *const c_char;
    fn proj_is_equivalent_to(left: *const c_void, right: *const c_void, criterion: c_int) -> c_int;
    fn proj_as_proj_string(
        ctx: *mut c_void,
        p: *const c_void,
        format: c_int,
        options: *const *const c_char,
    ) -> *const c_char;
    fn proj_coordoperation_get_accuracy(ctx: *mut c_void, p: *const c_void) -> f64;
    fn proj_get_area_of_use(
        ctx: *mut c_void,
        p: *const c_void,
        west: *mut f64,
        south: *mut f64,
        east: *mut f64,
        north: *mut f64,
        name: *mut *const c_char,
    ) -> c_int;
    fn proj_coordoperation_get_grid_used_count(ctx: *mut c_void, p: *const c_void) -> c_int;
    fn proj_coordoperation_get_grid_used(
        ctx: *mut c_void,
        p: *const c_void,
        index: c_int,
        short: *mut *const c_char,
        full: *mut *const c_char,
        package: *mut *const c_char,
        url: *mut *const c_char,
        direct: *mut c_int,
        open: *mut c_int,
        available: *mut c_int,
    ) -> c_int;
}
unsafe extern "C" fn destroy_operation(p: *mut c_void) {
    // SAFETY: Handle transfers the sole live PROJ allocation here.
    unsafe {
        proj_destroy(p);
    }
}
unsafe extern "C" fn destroy_context(p: *mut c_void) {
    // SAFETY: All dependent operations have already been dropped.
    unsafe {
        proj_context_destroy(p);
    }
}
fn failure(ctx: &Handle) -> Error {
    // SAFETY: Context is live; error string is copied before return.
    unsafe {
        backend_error(
            "PROJ",
            copy_string(proj_context_errno_string(
                ctx.ptr(),
                proj_context_errno(ctx.ptr()),
            )),
        )
    }
}
pub fn version() -> String {
    // SAFETY: Static native version string is copied.
    unsafe { copy_string(proj_info().version) }
}
/// CRS equivalence ignoring geographic latitude/longitude axis ORDER only.
/// Does not treat different datums, projected units or vertical systems as aliases.
pub fn equivalent_crs(left: &str, right: &str) -> Result<bool> {
    let left = cstring(left)?;
    let right = cstring(right)?;
    // SAFETY: Context owns both CRS objects; all objects and native strings are live.
    unsafe {
        let context = Handle::new(
            proj_context_create(),
            destroy_context,
            "PROJ",
            "cannot allocate context",
        )?;
        proj_context_set_enable_network(context.ptr(), 0);
        let left = Handle::new(
            proj_create(context.ptr(), left.as_ptr()),
            destroy_operation,
            "PROJ",
            "invalid CRS",
        )
        .map_err(|_| failure(&context))?;
        let right = Handle::new(
            proj_create(context.ptr(), right.as_ptr()),
            destroy_operation,
            "PROJ",
            "invalid CRS",
        )
        .map_err(|_| failure(&context))?;
        Ok(proj_is_equivalent_to(left.ptr(), right.ptr(), 2) != 0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AreaOfInterest {
    bounds: [f64; 4],
}
impl AreaOfInterest {
    /// west > east represents an antimeridian-crossing box.
    pub fn new(west: f64, south: f64, east: f64, north: f64) -> Result<Self> {
        for v in [west, south, east, north] {
            finite(v, "area bound")?;
        }
        if !(-180.0..=180.0).contains(&west)
            || !(-180.0..=180.0).contains(&east)
            || !(-90.0..=90.0).contains(&south)
            || !(-90.0..=90.0).contains(&north)
            || south >= north
            || west == east
        {
            return Err(Error::InvalidInput("invalid geographic area".into()));
        }
        Ok(Self {
            bounds: [west, south, east, north],
        })
    }
    pub fn bounds(self) -> [f64; 4] {
        self.bounds
    }
}
#[derive(Debug, Clone, Default)]
pub struct TransformOptions {
    pub area: Option<AreaOfInterest>,
    /// Desired published operation accuracy, not a guarantee or output covariance.
    pub maximum_accuracy_m: Option<f64>,
    /// Overrides native search paths. Include the directory containing proj.db.
    pub search_paths: Vec<std::path::PathBuf>,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coordinate4 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub decimal_year: f64,
}
#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Forward,
    Inverse,
}
#[derive(Debug, Clone)]
pub struct GridResource {
    pub name: String,
    pub full_name: String,
    pub package: String,
    pub url: String,
    pub available: bool,
    pub open_license: bool,
}
#[derive(Debug, Clone)]
pub struct OperationReport {
    pub name: String,
    pub definition: String,
    pub accuracy_m: Option<f64>,
    pub area: Option<[f64; 4]>,
    pub area_name: String,
    pub grids: Vec<GridResource>,
}
pub struct TransformedCoordinate {
    pub coordinate: Coordinate4,
    pub operation: OperationReport,
}
pub struct ProjTransform {
    // Field order ensures operations are destroyed before their owning context.
    operation: Handle,
    context: Handle,
    maximum_accuracy_m: Option<f64>,
}
impl ProjTransform {
    pub fn new(source_crs: &str, target_crs: &str, options: TransformOptions) -> Result<Self> {
        let source = cstring(source_crs)?;
        let target = cstring(target_crs)?;
        if let Some(accuracy) = options.maximum_accuracy_m {
            finite(accuracy, "accuracy")?;
            if accuracy < 0.0 {
                return Err(Error::InvalidInput("accuracy must be nonnegative".into()));
            }
        }
        let paths = options
            .search_paths
            .iter()
            .map(|p| {
                p.to_str()
                    .ok_or_else(|| Error::InvalidInput("PROJ path must be UTF-8".into()))
                    .and_then(cstring)
            })
            .collect::<Result<Vec<_>>>()?;
        let path_count = c_int::try_from(paths.len())
            .map_err(|_| Error::InvalidInput("too many search paths".into()))?;
        // SAFETY: Native ABI matches proj.h; all strings/option arrays live through calls.
        // All owned objects get the corresponding RAII destructor on every error path.
        unsafe {
            let info = proj_info();
            if (info.major, info.minor) < (9, 2) {
                return Err(backend_error("PROJ", "PROJ >= 9.2 required"));
            }
            let context = Handle::new(
                proj_context_create(),
                destroy_context,
                "PROJ",
                "cannot allocate context",
            )?;
            proj_context_set_enable_network(context.ptr(), 0);
            if !paths.is_empty() {
                let pointers: Vec<_> = paths.iter().map(|p| p.as_ptr()).collect();
                proj_context_set_search_paths(context.ptr(), path_count, pointers.as_ptr());
            }
            let source = Handle::new(
                proj_create(context.ptr(), source.as_ptr()),
                destroy_operation,
                "PROJ",
                "invalid source CRS/database unavailable",
            )
            .map_err(|_| failure(&context))?;
            let target = Handle::new(
                proj_create(context.ptr(), target.as_ptr()),
                destroy_operation,
                "PROJ",
                "invalid target CRS/database unavailable",
            )
            .map_err(|_| failure(&context))?;
            let area = options
                .area
                .map(|bounds| -> Result<Handle> {
                    let area = Handle::new(
                        proj_area_create(),
                        proj_area_destroy,
                        "PROJ",
                        "cannot allocate area",
                    )?;
                    let [w, s, e, n] = bounds.bounds;
                    proj_area_set_bbox(area.ptr(), w, s, e, n);
                    Ok(area)
                })
                .transpose()?;
            let mut option_strings = vec![cstring("ALLOW_BALLPARK=NO")?, cstring("ONLY_BEST=YES")?];
            if let Some(accuracy) = options.maximum_accuracy_m {
                option_strings.push(cstring(&format!("ACCURACY={accuracy}"))?);
            }
            let mut pointers: Vec<_> = option_strings.iter().map(|s| s.as_ptr()).collect();
            pointers.push(ptr::null());
            let operation = Handle::new(
                proj_create_crs_to_crs_from_pj(
                    context.ptr(),
                    source.ptr(),
                    target.ptr(),
                    area.as_ref().map_or(ptr::null_mut(), Handle::ptr),
                    pointers.as_ptr(),
                ),
                destroy_operation,
                "PROJ",
                "cannot create strict operation",
            )
            .map_err(|_| failure(&context))?;
            let normalized = Handle::new(
                proj_normalize_for_visualization(context.ptr(), operation.ptr()),
                destroy_operation,
                "PROJ",
                "cannot normalize axes",
            )
            .map_err(|_| failure(&context))?;
            // Locals (including the unnormalized operation) drop before context is moved.
            Ok(Self {
                operation: normalized,
                context,
                maximum_accuracy_m: options.maximum_accuracy_m,
            })
        }
    }
    /// Returns the operation actually used for THIS coordinate, not a generic CRS label.
    /// Accuracy/area are metadata; uncertainty and area enforcement remain caller work.
    pub fn transform(
        &mut self,
        coordinate: Coordinate4,
        direction: Direction,
    ) -> Result<TransformedCoordinate> {
        for value in [
            coordinate.x,
            coordinate.y,
            coordinate.z,
            coordinate.decimal_year,
        ] {
            finite(value, "coordinate")?;
        }
        // SAFETY: Owned context/operation are live; Coordinate matches PJ_COORD's ABI.
        // &mut self prevents concurrent use of its native context. Metadata is copied.
        unsafe {
            proj_errno_reset(self.operation.ptr());
            let output = proj_trans(
                self.operation.ptr(),
                match direction {
                    Direction::Forward => 1,
                    Direction::Inverse => -1,
                },
                Coordinate {
                    values: [
                        coordinate.x,
                        coordinate.y,
                        coordinate.z,
                        coordinate.decimal_year,
                    ],
                },
            )
            .values;
            let error = proj_errno(self.operation.ptr());
            if error != 0 {
                return Err(backend_error(
                    "PROJ",
                    copy_string(proj_context_errno_string(self.context.ptr(), error)),
                ));
            }
            for value in output {
                finite(value, "transformed coordinate")?;
            }
            let used = proj_trans_get_last_used_operation(self.operation.ptr());
            let used = if used.is_null() {
                None
            } else {
                Some(Handle::new(
                    used,
                    destroy_operation,
                    "PROJ",
                    "missing selected operation",
                )?)
            };
            let operation = self.report(used.as_ref().map_or(self.operation.ptr(), Handle::ptr));
            if let Some(maximum) = self.maximum_accuracy_m
                && operation
                    .accuracy_m
                    .is_none_or(|accuracy| accuracy > maximum)
            {
                return Err(backend_error(
                    "PROJ",
                    "selected operation has unknown or insufficient published accuracy",
                ));
            }
            Ok(TransformedCoordinate {
                coordinate: Coordinate4 {
                    x: output[0],
                    y: output[1],
                    z: output[2],
                    decimal_year: output[3],
                },
                operation,
            })
        }
    }
    unsafe fn report(&self, operation: *mut c_void) -> OperationReport {
        // SAFETY: Caller supplies a live operation in this context; all output pointers
        // reference stack storage and all returned strings are copied immediately.
        unsafe {
            let ctx = self.context.ptr();
            let accuracy = proj_coordoperation_get_accuracy(ctx, operation);
            let mut bounds = [0.0; 4];
            let mut name = ptr::null();
            let area_found = proj_get_area_of_use(
                ctx,
                operation,
                &mut bounds[0],
                &mut bounds[1],
                &mut bounds[2],
                &mut bounds[3],
                &mut name,
            ) != 0;
            let mut grids = Vec::new();
            for index in 0..proj_coordoperation_get_grid_used_count(ctx, operation) {
                let (mut short, mut full, mut package, mut url) =
                    (ptr::null(), ptr::null(), ptr::null(), ptr::null());
                let (mut direct, mut open, mut available) = (0, 0, 0);
                if proj_coordoperation_get_grid_used(
                    ctx,
                    operation,
                    index,
                    &mut short,
                    &mut full,
                    &mut package,
                    &mut url,
                    &mut direct,
                    &mut open,
                    &mut available,
                ) != 0
                {
                    grids.push(GridResource {
                        name: copy_string(short),
                        full_name: copy_string(full),
                        package: copy_string(package),
                        url: copy_string(url),
                        available: available != 0,
                        open_license: open != 0,
                    });
                }
            }
            OperationReport {
                name: copy_string(proj_get_name(operation)),
                definition: copy_string(proj_as_proj_string(ctx, operation, 0, ptr::null())),
                accuracy_m: (accuracy.is_finite() && accuracy >= 0.0).then_some(accuracy),
                area: area_found.then_some(bounds),
                area_name: copy_string(name),
                grids,
            }
        }
    }
}
