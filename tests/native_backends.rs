#![cfg(any(feature = "proj-backend", feature = "gdal-backend"))]
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "geoid-native-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[cfg(feature = "proj-backend")]
mod projection {
    use super::*;
    use geoid_atlas::backends::proj::*;
    fn point(x: f64, y: f64, z: f64) -> Coordinate4 {
        Coordinate4 {
            x,
            y,
            z,
            decimal_year: 2026.0,
        }
    }
    #[test]
    fn epsg_utm_axes_control_and_inverse_report() {
        let mut transform =
            ProjTransform::new("EPSG:4326", "EPSG:32613", TransformOptions::default()).unwrap();
        let result = transform
            .transform(point(-105.0, 40.0, 123.0), Direction::Forward)
            .unwrap();
        close(result.coordinate.x, 500_000.0, 1e-6);
        close(result.coordinate.y, 4_427_757.218_738, 0.001);
        close(result.coordinate.z, 123.0, 0.0); // 2D CRS passes z; not a vertical conversion
        assert!(!result.operation.definition.is_empty());
        assert!(!result.operation.name.is_empty());
        let round = transform
            .transform(result.coordinate, Direction::Inverse)
            .unwrap();
        close(round.coordinate.x, -105.0, 1e-10);
        close(round.coordinate.y, 40.0, 1e-10);
        close(round.coordinate.decimal_year, 2026.0, 0.0);
    }
    #[test]
    fn authoritative_unit_conversion_and_3d_geocentric_height() {
        let mut units =
            ProjTransform::new("EPSG:2230", "EPSG:26946", TransformOptions::default()).unwrap();
        let result = units
            .transform(
                point(4_760_096.421_921, 3_744_293.729_449, 0.0),
                Direction::Forward,
            )
            .unwrap();
        close(result.coordinate.x, 1_450_880.291_060_502_2, 1e-7);
        close(result.coordinate.y, 1_141_263.011_160_478_2, 1e-7);
        let mut cartesian =
            ProjTransform::new("EPSG:4979", "EPSG:4978", TransformOptions::default()).unwrap();
        let result = cartesian
            .transform(point(0.0, 0.0, 100.0), Direction::Forward)
            .unwrap();
        close(result.coordinate.x, 6_378_237.0, 1e-8);
        close(result.coordinate.y, 0.0, 1e-8);
        close(result.coordinate.z, 0.0, 1e-8);
    }
    #[test]
    fn malformed_crs_inputs_accuracy_and_projection_domain_fail() {
        assert!(
            ProjTransform::new(
                "EPSG:does-not-exist",
                "EPSG:4326",
                TransformOptions::default()
            )
            .is_err()
        );
        assert!(
            ProjTransform::new("EPSG:4326\0", "EPSG:3857", TransformOptions::default()).is_err()
        );
        assert!(
            ProjTransform::new(
                "EPSG:4326",
                "EPSG:3857",
                TransformOptions {
                    maximum_accuracy_m: Some(-1.0),
                    ..Default::default()
                }
            )
            .is_err()
        );
        let mut transform =
            ProjTransform::new("EPSG:4326", "EPSG:3857", TransformOptions::default()).unwrap();
        assert!(
            transform
                .transform(point(0.0, 95.0, 0.0), Direction::Forward)
                .is_err()
        );
        assert!(
            transform
                .transform(point(f64::NAN, 40.0, 0.0), Direction::Forward)
                .is_err()
        );
        assert!(AreaOfInterest::new(-180.0, 20.0, 180.0, 10.0).is_err());
        assert!(AreaOfInterest::new(170.0, -20.0, -170.0, 20.0).is_ok());
        // Failed calls must not poison later coordinates.
        assert!(
            transform
                .transform(point(0.0, 0.0, 0.0), Direction::Forward)
                .is_ok()
        );
    }
    #[test]
    fn dynamic_itrf_operation_uses_coordinate_epoch_and_round_trips() {
        let mut transform =
            ProjTransform::new("EPSG:5332", "EPSG:7789", TransformOptions::default()).unwrap();
        let initial = point(3_657_660.66, 255_768.55, 5_201_382.11);
        let at_2010 = transform
            .transform(
                Coordinate4 {
                    decimal_year: 2010.0,
                    ..initial
                },
                Direction::Forward,
            )
            .unwrap();
        let at_2020 = transform
            .transform(
                Coordinate4 {
                    decimal_year: 2020.0,
                    ..initial
                },
                Direction::Forward,
            )
            .unwrap();
        assert!(at_2020.operation.definition.contains("helmert"));
        assert!(at_2020.operation.definition.contains("t_epoch"));
        let displacement = ((at_2020.coordinate.x - at_2010.coordinate.x).powi(2)
            + (at_2020.coordinate.y - at_2010.coordinate.y).powi(2)
            + (at_2020.coordinate.z - at_2010.coordinate.z).powi(2))
        .sqrt();
        assert!(displacement > 0.0001 && displacement < 0.02);
        let round = transform
            .transform(at_2020.coordinate, Direction::Inverse)
            .unwrap();
        close(round.coordinate.x, initial.x, 1e-8);
        close(round.coordinate.y, initial.y, 1e-8);
        close(round.coordinate.z, initial.z, 1e-8);
        close(round.coordinate.decimal_year, 2020.0, 0.0);
    }
    #[test]
    fn missing_vertical_grid_is_an_error_in_database_only_context() {
        let fixture = Fixture::new();
        let data = std::env::var_os("GEOID_TEST_PROJ_DATA_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("/usr/share/proj"));
        std::fs::copy(data.join("proj.db"), fixture.0.join("proj.db"))
            .expect("set GEOID_TEST_PROJ_DATA_DIR to installed proj.db directory");
        let transform = ProjTransform::new(
            "EPSG:4326+5773",
            "EPSG:4979",
            TransformOptions {
                search_paths: vec![fixture.0.clone()],
                ..Default::default()
            },
        );
        match transform {
            Err(_) => (),
            Ok(mut transform) => assert!(
                transform
                    .transform(point(-105.0, 40.0, 100.0), Direction::Forward)
                    .is_err(),
                "EGM96 conversion must not silently pass through height"
            ),
        }
    }
    #[test]
    fn area_and_accuracy_filter_forbid_ballpark_operations() {
        let options = TransformOptions {
            area: Some(AreaOfInterest::new(-110.0, 35.0, -100.0, 45.0).unwrap()),
            maximum_accuracy_m: Some(0.0),
            ..Default::default()
        };
        // A genuine NAD27/NAD83 datum shift is not a zero-error conversion.
        match ProjTransform::new("EPSG:4267", "EPSG:4269", options) {
            Err(_) => (),
            Ok(mut transform) => assert!(
                transform
                    .transform(point(-105.0, 40.0, 0.0), Direction::Forward)
                    .is_err()
            ),
        }
    }
}

#[cfg(feature = "gdal-backend")]
mod imagery {
    use super::*;
    use geoid_atlas::Error;
    use geoid_atlas::backends::gdal::*;
    fn fixture() -> Fixture {
        let fixture = Fixture::new();
        std::fs::write(fixture.0.join("source.asc"),"ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n1 -9999 3\n4 5 6\n").unwrap();
        std::fs::write(
            fixture.0.join("mask.asc"),
            "ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\n255 255 0\n255 255 255\n",
        )
        .unwrap();
        std::fs::write(fixture.0.join("source.vrt"),r#"<VRTDataset rasterXSize="3" rasterYSize="2">
<SRS>EPSG:4326</SRS><GeoTransform>-105,0.01,0.002,40,0.001,-0.01</GeoTransform>
<VRTRasterBand dataType="Float64" band="1"><NoDataValue>-9999</NoDataValue><Scale>2</Scale><Offset>10</Offset><UnitType>m</UnitType>
<SimpleSource><SourceFilename relativeToVRT="1">source.asc</SourceFilename><SourceBand>1</SourceBand><SrcRect xOff="0" yOff="0" xSize="3" ySize="2"/><DstRect xOff="0" yOff="0" xSize="3" ySize="2"/></SimpleSource>
<MaskBand><VRTRasterBand dataType="Byte"><SimpleSource><SourceFilename relativeToVRT="1">mask.asc</SourceFilename><SourceBand>1</SourceBand><SrcRect xOff="0" yOff="0" xSize="3" ySize="2"/><DstRect xOff="0" yOff="0" xSize="3" ySize="2"/></SimpleSource></VRTRasterBand></MaskBand>
</VRTRasterBand></VRTDataset>"#).unwrap();
        fixture
    }
    #[test]
    fn reads_rotated_window_crs_mask_nodata_and_physical_units() {
        let fixture = fixture();
        let mut raster = GdalRaster::open(fixture.0.join("source.vrt")).unwrap();
        assert_eq!(raster.metadata().width, 3);
        assert!(raster.metadata().crs_wkt.contains("WGS 84"));
        let window = raster.read_window(1, 0, 0, 3, 2).unwrap();
        assert_eq!(window.unit, "m");
        assert_eq!(
            window.cells,
            vec![Some(12.0), None, None, Some(18.0), Some(20.0), Some(22.0)]
        );
        let (x, y) = window.transform.world(0.0, 0.0);
        close(x, -104.994, 1e-12);
        close(y, 39.9955, 1e-12);
        assert_eq!(window.nearest(x, y), Ok(12.0));
        let subwindow = raster.read_window(1, 1, 1, 2, 1).unwrap();
        assert_eq!(subwindow.cells, vec![Some(20.0), Some(22.0)]);
        let expected = window.transform.world(1.0, 1.0);
        assert_eq!(subwindow.transform.world(0.0, 0.0), expected);
        let missing = window.transform.world(1.0, 0.0);
        assert_eq!(window.nearest(missing.0, missing.1), Err(Error::NoData));
        assert!(raster.read_window(0, 0, 0, 1, 1).is_err());
        assert!(raster.read_window(1, usize::MAX, 0, 1, 1).is_err());
        assert!(GdalRaster::open(fixture.0.join("source.asc")).is_err()); // absent CRS
    }
    #[test]
    fn warps_real_geotiff_and_refuses_existing_destination() {
        let fixture = fixture();
        let output = fixture.0.join("warped");
        let merc_x = |lon: f64| 6_378_137.0 * lon.to_radians();
        let merc_y = |lat: f64| {
            6_378_137.0
                * (std::f64::consts::FRAC_PI_4 + lat.to_radians() / 2.0)
                    .tan()
                    .ln()
        };
        let request = || WarpRequest {
            target_crs: "EPSG:3857",
            bounds: [
                merc_x(-105.0),
                merc_y(39.98),
                merc_x(-104.966),
                merc_y(40.003),
            ],
            width: 8,
            height: 6,
            resampling: Resampling::Nearest,
            source_decimal_year: None,
            target_decimal_year: None,
        };
        let path = warp_to_new_directory(fixture.0.join("source.vrt"), &output, request()).unwrap();
        let mut raster = GdalRaster::open(&path).unwrap();
        assert_eq!(raster.metadata().width, 8);
        assert_eq!(raster.metadata().height, 6);
        assert!(raster.metadata().crs_wkt.contains("3857"));
        let cells = raster.read_window(1, 0, 0, 8, 6).unwrap().cells;
        assert!(cells.iter().any(Option::is_some));
        assert!(cells.iter().any(Option::is_none));
        let bytes = std::fs::read(&path).unwrap();
        assert!(warp_to_new_directory(fixture.0.join("source.vrt"), &output, request()).is_err());
        assert_eq!(std::fs::read(path).unwrap(), bytes);
        #[cfg(feature = "proj-backend")]
        {
            use geoid_atlas::coordinates::Crs;
            let window = raster.read_window(1, 0, 0, 8, 6).unwrap();
            assert!(window.clone().into_core_raster(Crs::Wgs84).is_err());
            assert!(window.into_core_raster(Crs::WebMercator).is_ok());
        }
    }
}
