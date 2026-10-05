use geoid_atlas::Result;
use geoid_atlas::backends::{gdal, proj};
use proj::{Coordinate4, Direction, ProjTransform, TransformOptions};

fn main() -> Result<()> {
    println!("Native PROJ {}; GDAL {}", proj::version(), gdal::version());
    let mut transform = ProjTransform::new("EPSG:4326", "EPSG:32613", TransformOptions::default())?;
    let result = transform.transform(
        Coordinate4 {
            x: -105.0,
            y: 40.0,
            z: 123.0,
            decimal_year: 2026.0,
        },
        Direction::Forward,
    )?;
    println!(
        "WGS84 -> UTM13N: {:.3}, {:.3}; z {:.3}",
        result.coordinate.x, result.coordinate.y, result.coordinate.z
    );
    println!("Selected operation: {}", result.operation.name);
    println!("Pipeline: {}", result.operation.definition);
    println!(
        "Published accuracy: {:?}; grids: {:?}",
        result.operation.accuracy_m, result.operation.grids
    );
    assert!((result.coordinate.x - 500_000.0).abs() < 1e-6);
    assert!((result.coordinate.y - 4_427_757.218_738).abs() < 0.001);
    // Optional read-only local GeoTIFF/VRT/etc. inspection; no implicit download/warp.
    if let Some(path) = std::env::args_os().nth(1) {
        let mut raster = gdal::GdalRaster::open(path)?;
        println!("Raster: {:?}", raster.metadata());
        let width = raster.metadata().width.min(4);
        let height = raster.metadata().height.min(4);
        println!(
            "First band window: {:?}",
            raster.read_window(1, 0, 0, width, height)?.cells
        );
    }
    Ok(())
}
