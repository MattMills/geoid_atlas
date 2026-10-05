# Authoritative native geospatial backends

The default Rust core remains dependency-free and builds without native GIS software. Two opt-in Cargo features link the public C APIs of **PROJ >= 9.2** and **GDAL >= 3.8**. There are no additional Rust crates. The adapters use owned native handles with matching destructors; each context/dataset is kept private and is not declared Send/Sync. Use one transformation per worker rather than sharing native contexts. The small FFI declarations follow the upstream public headers; installing development libraries supplies the linker files without generating Rust bindings.

## Install and enable

On Ubuntu 24.04 or another supported Debian-family system:

```sh
sudo apt-get install libproj-dev proj-data libgdal-dev
export PROJ_NETWORK=OFF
cargo test --locked --all-features --all-targets
cargo test --locked --all-features --doc
cargo run --locked --all-features --example native_geospatial
```

For macOS, install `proj` and `gdal` through the native package manager and expose their library directories to the linker. Windows needs compatible native libraries/import libraries and runtime DLLs. Those platforms are not validated by this change; Linux native integration is exercised in CI. `GEOID_NATIVE_LIB_DIR` optionally adds a nonstandard native link directory at build time. Runtime dynamic-library search paths are configured separately by the host. A native runtime version check cannot replace installing ABI-compatible libraries.

PROJ requires its database (`proj.db`) and any operation's grids. Install grid files separately with their provenance/licenses and checksums. `PROJ_DATA` selects the normal native data path. `TransformOptions::search_paths` overrides that for an individual PROJ context; include the database directory as well as grid directories. Tests that isolate missing grids copy only `proj.db`; `GEOID_TEST_PROJ_DATA_DIR` identifies its source directory (defaults to `/usr/share/proj`). Nothing in the adapter downloads a database or grid.

## PROJ coordinate transformations

```rust
use geoid_atlas::backends::proj::*;
let mut operation = ProjTransform::new(
    "EPSG:4326", "EPSG:32613", TransformOptions::default(),
)?;
let result = operation.transform(Coordinate4 {
    x: -105.0, y: 40.0, z: 123.0, decimal_year: 2026.0,
}, Direction::Forward)?;
assert!((result.coordinate.x - 500_000.0).abs() < 1e-6);
println!("{}", result.operation.definition);
# Ok::<(), geoid_atlas::Error>(())
```

Source/target definitions can be EPSG identifiers, compound CRS identifiers, WKT, or other CRS definitions accepted by the installed PROJ database. An arbitrary operation pipeline is not a CRS definition. Axis order is normalized for visualization: geographic coordinates are **longitude, latitude** in degrees; projected/geocentric coordinates use the actual CRS axes/units. Coordinates are not globally forced to metres, because some authoritative CRSs use survey feet. `Coordinate4::z` follows the declared CRS, and `decimal_year` is an explicit coordinate epoch in the operation's convention, not seconds from the library's `Epoch`. These APIs transform coordinate events, not velocity vectors or state covariance. Use the core frame/kinematic APIs or a suitable Jacobian and time derivative for those.

Use **3D or compound** source/target CRSs when transforming heights, such as EPSG:4979 versus EPSG:4978, or `EPSG:4326+5773` for WGS84 with EGM96 orthometric height. A 2D CRS pair can pass z unchanged; supplying a number in z does not make that a vertical transformation. A coordinate epoch is not automatic physical propagation of a station to another observation time. Coordinate-metadata objects specifying distinct source/target epochs and proper-time conversion are not exposed by this adapter yet.

All transformations use `ALLOW_BALLPARK=NO`, `ONLY_BEST=YES`, and context-local grid networking disabled. Missing best-operation grids cause a native error during creation or execution, rather than silently using a ballpark identity. Optional desired accuracy filters candidate operations and rejects outputs whose actual published accuracy is unknown or exceeds the requested maximum. A geographic area of interest helps select an operation; it does not clip coordinates or prove accuracy outside its area of use.

`TransformedCoordinate::operation` copies the operation actually selected for that point: name, pipeline, published accuracy where known, area of use, and referenced grid names/paths/URLs/availability/license flags. Preserve this along with the PROJ version and database/grid versions for reproducibility. PROJ's published accuracy is not measurement noise, a full covariance, or a guaranteed error bound. Different points can select different operations. `equivalent_crs` checks native CRS equivalence while ignoring geographic axis-order differences; it does not equate distinct datums or units.

## GDAL raster ingestion and warping

`GdalRaster::open` reads an existing local raster file supported by the installed GDAL drivers (including GeoTIFF and VRT). It requires an explicit CRS and affine georeferencing. It does not synthesize georeferencing from a raw photograph, decode camera attitude, or ingest an arbitrary raster without a CRS. GDAL formats can reference sidecar/source files; retain those dependencies with the data.

`read_window` reads a bounded, 1-based band into a `RasterWindow` with its WKT, pixel-centre affine transform and native band unit label (which may be unspecified). It preserves masked/nodata cells as None and applies band scale/offset to finite values. Window coordinates are in the original CRS. Masks are treated as validity (nonzero means valid), not fractional coverage weights. The adapter reads scalar band values, not display-ready RGB color management or an image viewer. No vertical-datum meaning or metre units are inferred from a DEM file's numeric band. Convert values and uncertainties to the core quantity's SI units explicitly before registering a field in the atlas.

With both features enabled, `window.into_core_raster(crs)` attaches an already-compatible WGS84 geographic, Web Mercator, or WGS84 UTM window to the existing raster/atlas APIs **only after native CRS equivalence validation**. A mismatch returns an error. Warp first if needed; this method cannot relabel a NAD83 raster as WGS84. Planetary rasters need an explicit body-specific adapter. A DEM can then wrap the resulting raster with its explicit `VerticalDatum` and enter the source-backed atlas.

```rust
use geoid_atlas::backends::gdal::*;
// Existing georeferenced input; target bounds and dimensions are explicit.
# fn demonstrate(input: &std::path::Path, new_output_dir: &std::path::Path) -> geoid_atlas::Result<()> {
let output = warp_to_new_directory(input, new_output_dir, WarpRequest {
    target_crs: "EPSG:3857",
    bounds: [-11_690_000.0, 4_860_000.0, -11_680_000.0, 4_870_000.0],
    width: 1024, height: 1024,
    resampling: Resampling::Bilinear,
    source_decimal_year: None, target_decimal_year: None,
})?;
let mut image = GdalRaster::open(output)?;
let window = image.read_window(1, 0, 0, 256, 256)?;
assert_eq!(window.cells.len(), 256 * 256);
# Ok(()) }
```

The warp writes `raster.tif` inside a **new directory**, refusing an existing directory to avoid updating user imagery. A failed warp can leave partial output there; it is not deleted automatically. It uses GDAL's warper with explicit resampling, target bounds/dimensions, `-et 0` exact transformer evaluation, and strict PROJ operation selection. This removes the transformer's approximation threshold, but does not guarantee image resolution, interpolation accuracy or a fully precise datum model. The installed GDAL must link PROJ >= 9.2. GDAL's PROJ network setting is application-global: the adapter rejects enabled networking and does not change another application's global setting. Launch with `PROJ_NETWORK=OFF` and do not toggle native GIS configuration concurrently with operations.

Optional source/target decimal coordinate epochs are passed through GDAL's native flags; combinations depend on the linked PROJ version and supported velocity models/grids. No raw-camera orthorectification, automatic vertical conversion policy, uncertainty resampling, antimeridian splitting, or arbitrary pipeline override is inferred. Warping an elevation band may invoke GDAL's vertical shift behavior when a compound CRS supports it; choose the declared height reference and inspect the resulting units. A raster warp may choose different operations across pixels; the point-transform report is not a provenance report for an entire warped image. Record source/target definitions, epochs, grid installations, native versions and warp settings.

## Licensing

PROJ is permissively **MIT-style licensed**, allowing modification, commercial use and redistribution while retaining its copyright/permission notices. GDAL's general license is likewise **MIT-style**, with additional notices for individual components. This is compatible with keeping `geoid-atlas` under MIT; enabling the backends does not change this crate's license.

Verified upstream sources:

- [PROJ 9.6.0 COPYING](https://github.com/OSGeo/PROJ/blob/9.6.0/COPYING), retained in `licenses/PROJ-9.6.0-COPYING.txt`.
- [GDAL 3.10.3 LICENSE.TXT](https://github.com/OSGeo/gdal/blob/v3.10.3/LICENSE.TXT), retained in `licenses/GDAL-3.10.3-LICENSE.txt`.
- [EPSG dataset terms](https://epsg.org/terms-of-use.html) and [PROJ-data](https://github.com/OSGeo/PROJ-data) for the database/grid resource provenance.

The repository does not vendor native binaries, EPSG databases, or correction grids. When bundling those, preserve the notices and terms for the versions and datasets actually distributed. GDAL optional drivers/native dependencies and individual grid datasets can carry additional licenses; select and review the native distribution you ship. A PROJ grid's `open_license` flag is metadata, not a substitute for its full license.
