use geoid_atlas::Error;
use geoid_atlas::atlas::*;
use geoid_atlas::coordinates::*;
use geoid_atlas::frames::*;
use geoid_atlas::raster::*;

fn epoch(t: f64) -> Epoch {
    Epoch::new(t).unwrap()
}
fn point() -> Geodetic {
    Geodetic::new(0.0, 0.0, -100.0).unwrap()
}
fn atlas() -> Atlas {
    let mut frames = FrameGraph::new("barycentric");
    let root = frames.root();
    let earth = frames.add("earth-fixed", root, Pose::IDENTITY).unwrap();
    let mars = frames
        .add(
            "mars-fixed",
            root,
            Pose {
                origin: State::stationary(Vec3::new(1e11, 0.0, 0.0).unwrap()),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let mut atlas = Atlas::new(frames);
    atlas
        .add_body(
            "earth",
            Body {
                name: "Earth".into(),
                fixed_frame: earth,
                ellipsoid: Ellipsoid::WGS84,
                gravitational_parameter_m3_s2: 3.986e14,
            },
        )
        .unwrap();
    atlas
        .add_body(
            "mars",
            Body {
                name: "Mars".into(),
                fixed_frame: mars,
                ellipsoid: Ellipsoid::new(3_396_190.0, 1.0 / 169.8).unwrap(),
                gravitational_parameter_m3_s2: 4.282e13,
            },
        )
        .unwrap();
    atlas
}
fn source(atlas: &mut Atlas, id: &str, group: &str, kind: SourceKind) {
    atlas
        .add_source(Source {
            id: id.into(),
            title: id.into(),
            reference: "synthetic test fixture".into(),
            kind,
            independence_group: group.into(),
        })
        .unwrap();
}
fn field(value: f64, std_dev: f64) -> ConstantField {
    ConstantField {
        estimate: Estimate::new(value, std_dev).unwrap(),
    }
}

#[test]
fn independent_fusion_keeps_correlated_products_from_double_counting() {
    let mut atlas = atlas();
    source(&mut atlas, "survey", "a", SourceKind::Observation);
    source(&mut atlas, "resampled", "a", SourceKind::Observation);
    source(&mut atlas, "other", "b", SourceKind::Observation);
    atlas
        .add_scalar(
            "earth",
            "survey",
            Quantity::EllipsoidalElevationM,
            None,
            field(100.0, 2.0),
        )
        .unwrap();
    atlas
        .add_scalar(
            "earth",
            "resampled",
            Quantity::EllipsoidalElevationM,
            None,
            field(100.0, 3.0),
        )
        .unwrap();
    atlas
        .add_scalar(
            "earth",
            "other",
            Quantity::EllipsoidalElevationM,
            None,
            field(100.0, 2.0),
        )
        .unwrap();
    let e = atlas
        .expectation(
            "earth",
            Quantity::EllipsoidalElevationM,
            point(),
            epoch(0.0),
            EvidenceSelection::Observed,
        )
        .unwrap();
    assert_eq!(e.estimate.value(), 100.0);
    assert!((e.estimate.std_dev() - 2.0_f64.sqrt()).abs() < 1e-12);
    assert_eq!(e.contributing_sources, vec!["survey", "other"]);
    assert_eq!(e.redundant_sources, vec!["resampled"]);
}

#[test]
fn body_time_and_observation_scopes_are_enforced() {
    let mut atlas = atlas();
    source(&mut atlas, "observation", "a", SourceKind::Observation);
    source(&mut atlas, "model", "b", SourceKind::Model);
    atlas
        .add_scalar(
            "mars",
            "observation",
            Quantity::ConductivitySPerM,
            Some(Validity::new(epoch(0.0), epoch(10.0)).unwrap()),
            field(1.0, 0.1),
        )
        .unwrap();
    atlas
        .add_scalar(
            "mars",
            "model",
            Quantity::ConductivitySPerM,
            None,
            field(9.0, 0.1),
        )
        .unwrap();
    let e = atlas
        .expectation(
            "mars",
            Quantity::ConductivitySPerM,
            point(),
            epoch(5.0),
            EvidenceSelection::Observed,
        )
        .unwrap();
    assert_eq!(e.estimate.value(), 1.0);
    assert!(matches!(
        atlas.expectation(
            "earth",
            Quantity::ConductivitySPerM,
            point(),
            epoch(5.0),
            EvidenceSelection::All
        ),
        Err(Error::NoData)
    ));
    assert!(matches!(
        atlas.expectation(
            "mars",
            Quantity::ConductivitySPerM,
            point(),
            epoch(11.0),
            EvidenceSelection::Observed
        ),
        Err(Error::NoData)
    ));
    assert_eq!(
        atlas
            .expectation(
                "mars",
                Quantity::ConductivitySPerM,
                point(),
                epoch(11.0),
                EvidenceSelection::Modeled
            )
            .unwrap()
            .estimate
            .value(),
        9.0
    );
}

#[test]
fn disagreement_inflates_uncertainty_and_unknown_sources_fail() {
    let mut atlas = atlas();
    source(&mut atlas, "a", "a", SourceKind::Observation);
    source(&mut atlas, "b", "b", SourceKind::Observation);
    atlas
        .add_scalar(
            "earth",
            "a",
            Quantity::EllipsoidalElevationM,
            None,
            field(0.0, 1.0),
        )
        .unwrap();
    atlas
        .add_scalar(
            "earth",
            "b",
            Quantity::EllipsoidalElevationM,
            None,
            field(10.0, 1.0),
        )
        .unwrap();
    let e = atlas
        .expectation(
            "earth",
            Quantity::EllipsoidalElevationM,
            point(),
            epoch(0.0),
            EvidenceSelection::All,
        )
        .unwrap();
    assert_eq!(e.estimate.value(), 5.0);
    assert!((e.estimate.std_dev() - 5.0).abs() < 1e-12);
    assert_eq!(e.reduced_chi_squared, Some(50.0));
    assert!(matches!(
        atlas.add_scalar(
            "earth",
            "missing",
            Quantity::GeoidUndulationM,
            None,
            field(0.0, 1.0)
        ),
        Err(Error::UnknownSource(_))
    ));
    assert!(
        atlas
            .add_scalar(
                "earth",
                "a",
                Quantity::EllipsoidalElevationM,
                None,
                field(1.0, 1.0)
            )
            .is_err()
    );
}

#[test]
fn vertical_conversion_site_and_planetary_raster_crs() {
    let mut atlas = atlas();
    source(&mut atlas, "dem", "a", SourceKind::Observation);
    source(&mut atlas, "geoid", "b", SourceKind::Model);
    atlas
        .add_scalar(
            "earth",
            "dem",
            Quantity::OrthometricElevationM,
            None,
            field(-100.0, 1.0),
        )
        .unwrap();
    atlas
        .add_scalar(
            "earth",
            "geoid",
            Quantity::GeoidUndulationM,
            None,
            field(30.0, 1.0),
        )
        .unwrap();
    let state = atlas
        .terrain_site_in_root(
            "earth",
            point(),
            epoch(0.0),
            VerticalDatum::Orthometric,
            EvidenceSelection::All,
        )
        .unwrap();
    assert!((state.position_m.x - (WGS84_A - 70.0)).abs() < 1e-9);
    let raster = Raster::new(
        1,
        1,
        Affine::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0).unwrap(),
        Crs::Wgs84,
        vec![Some(0.0)],
    )
    .unwrap();
    assert!(
        atlas
            .add_raster(
                "mars",
                "dem",
                Quantity::EllipsoidalElevationM,
                None,
                RasterField {
                    raster,
                    std_dev: 1.0
                }
            )
            .is_err()
    );
}

#[test]
fn drg_sampling_and_rf_bands_do_not_mix() {
    let mut atlas = atlas();
    source(&mut atlas, "map", "a", SourceKind::Observation);
    let rgb = Rgb {
        red: 10,
        green: 20,
        blue: 30,
    };
    let image = Raster::new(
        1,
        1,
        Affine::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0).unwrap(),
        Crs::Planetographic,
        vec![Some(rgb)],
    )
    .unwrap();
    atlas.add_drg("mars", "map", None, image).unwrap();
    assert_eq!(
        atlas.drg_pixel("mars", "map", point(), epoch(0.0)).unwrap(),
        rgb
    );
    let band = Band::new(1e9, 1.1e9).unwrap();
    atlas
        .add_scalar(
            "mars",
            "map",
            Quantity::ReceivedPowerDbm(band),
            None,
            field(-100.0, 3.0),
        )
        .unwrap();
    assert!(matches!(
        atlas.expectation(
            "mars",
            Quantity::ReceivedPowerDbm(Band::new(2e9, 2.1e9).unwrap()),
            point(),
            epoch(0.0),
            EvidenceSelection::All
        ),
        Err(Error::NoData)
    ));
}
