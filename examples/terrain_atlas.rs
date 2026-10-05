use geoid_atlas::Result;
use geoid_atlas::atlas::*;
use geoid_atlas::coordinates::*;
use geoid_atlas::frames::*;
use geoid_atlas::raster::*;

fn main() -> Result<()> {
    let mut frames = FrameGraph::new("Earth-centred inertial demo");
    let earth = frames.add("Earth-fixed demo", frames.root(), Pose::IDENTITY)?;
    let mut atlas = Atlas::new(frames);
    atlas.add_body(
        "earth",
        Body {
            name: "Earth".into(),
            fixed_frame: earth,
            ellipsoid: Ellipsoid::WGS84,
            gravitational_parameter_m3_s2: 3.986_004_418e14,
        },
    )?;
    for (id, kind) in [
        ("dem", SourceKind::Observation),
        ("geoid", SourceKind::Model),
    ] {
        atlas.add_source(Source {
            id: id.into(),
            title: format!("Synthetic {id}"),
            reference: "built-in example fixture; not real survey data".into(),
            kind,
            independence_group: id.into(),
        })?;
    }
    let dem = Dem::from_esri_ascii(
        "ncols 2\nnrows 2\nxllcenter -105.0\nyllcenter 40.0\ncellsize 0.01\n1600 1610\n1620 1630\n",
        Crs::Wgs84,
        VerticalDatum::Orthometric,
    )?;
    atlas.add_raster(
        "earth",
        "dem",
        Quantity::OrthometricElevationM,
        None,
        RasterField {
            raster: dem.raster,
            std_dev: 2.0,
        },
    )?;
    atlas.add_scalar(
        "earth",
        "geoid",
        Quantity::GeoidUndulationM,
        None,
        ConstantField {
            estimate: Estimate::new(-20.0, 0.5)?,
        },
    )?;
    let point = Geodetic::new(40.005, -104.995, 0.0)?;
    let epoch = Epoch::new(0.0)?;
    let height = atlas.expectation(
        "earth",
        Quantity::OrthometricElevationM,
        point,
        epoch,
        EvidenceSelection::Observed,
    )?;
    let site = atlas.terrain_site_in_root(
        "earth",
        point,
        epoch,
        VerticalDatum::Orthometric,
        EvidenceSelection::All,
    )?;
    assert!((height.estimate.value() - 1615.0).abs() < 1e-6);
    println!(
        "Orthometric elevation: {:.3} +/- {:.3} m",
        height.estimate.value(),
        height.estimate.std_dev()
    );
    println!("Root-frame site: {:?} m", site.position_m);
    Ok(())
}
