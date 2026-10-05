use geoid_atlas::coordinates::{Ecef, Ellipsoid, Geodetic, WGS84_A, WGS84_F};
use geoid_atlas::frames::Vec3;
use geoid_atlas::geodesics::{
    MEAN_EARTH_RADIUS_M, SurfaceGeodesic, great_circle_m, ground_distance_m,
};

fn point(lat: f64, lon: f64) -> Geodetic {
    Geodetic::from_validated(lat, lon, 0.0)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16} != {expected:.16}, error {}, tolerance {tolerance}",
        (actual - expected).abs()
    );
}
fn angle_error(actual: f64, expected: f64) -> f64 {
    ((actual - expected + 180.0).rem_euclid(360.0) - 180.0).abs()
}
fn rows(data: &str) -> Vec<Vec<f64>> {
    data.lines()
        .map(|line| {
            line.split_whitespace()
                .map(|v| v.parse().unwrap())
                .collect()
        })
        .collect()
}

#[test]
fn geographiclib_official_twenty_inverse_and_direct() {
    let geodesic = SurfaceGeodesic::wgs84();
    let cases = rows(include_str!("data/geodesics/GeographicLib-20.dat"));
    assert_eq!(cases.len(), 20);
    for v in cases {
        assert_eq!(v.len(), 12);
        let start = point(v[0], v[1]);
        let end = point(v[3], v[4]);
        let inv = geodesic.inverse(start, end);
        close(inv.distance_m, v[6], 1e-7);
        close(angle_error(inv.initial_azimuth_deg, v[2]), 0.0, 1e-12);
        close(angle_error(inv.final_azimuth_deg, v[5]), 0.0, 1e-12);
        close(inv.auxiliary_arc_deg, v[7], 5e-13);
        let direct = geodesic.direct(start, v[2], v[6]).unwrap();
        close(direct.point.latitude_deg(), v[3], 5e-13);
        close(angle_error(direct.point.longitude_deg(), v[4]), 0.0, 5e-13);
        close(angle_error(direct.final_azimuth_deg, v[5]), 0.0, 5e-13);
        // A negative distance follows the forward endpoint heading backwards.
        let reverse = geodesic.direct(end, v[5], -v[6]).unwrap();
        close(ground_distance_m(reverse.point, start), 0.0, 1e-7);
    }
}

#[test]
fn published_geodtest_hundred_inverse_and_direct() {
    let geodesic = SurfaceGeodesic::wgs84();
    let cases = rows(include_str!("data/geodesics/GeodTest-100.dat"));
    assert_eq!(cases.len(), 100);
    for (i, v) in cases.iter().enumerate() {
        assert_eq!(v.len(), 10);
        let start = point(v[0], v[1]);
        let end = point(v[3], v[4]);
        let inv = geodesic.inverse(start, end);
        close(inv.distance_m, v[6], 1e-7);
        close(inv.auxiliary_arc_deg, v[7], 2e-10);
        // Reduced length measures sensitivity of endpoint displacement to heading.
        // Rounded endpoints cannot resolve azimuth at a coincident/conjugate pair.
        for (actual, reference) in [
            (inv.initial_azimuth_deg, v[2]),
            (inv.final_azimuth_deg, v[5]),
        ] {
            let transverse_error = v[8].abs() * angle_error(actual, reference).to_radians();
            assert!(
                transverse_error < 2e-7,
                "row {} transverse error {transverse_error}",
                i + 1
            );
        }
        let direct = geodesic.direct(start, v[2], v[6]).unwrap();
        close(ground_distance_m(direct.point, end), 0.0, 1e-7);
        // Direct endpoint azimuth remains well-conditioned given the initial heading.
        close(angle_error(direct.final_azimuth_deg, v[5]), 0.0, 2e-8);
        let reverse = geodesic.direct(end, v[5], -v[6]).unwrap();
        close(ground_distance_m(reverse.point, start), 0.0, 1e-7);
    }
}

#[test]
fn published_antipodal_and_polar_regressions() {
    // GeographicLib GeodSolve6/9/10/11/99, preserving the published tolerances.
    for (lat1, lat2, lon2, distance) in [
        (
            88.202499451857,
            -88.202499451857,
            179.98102203299286,
            20003898.214,
        ),
        (
            89.262080389218,
            -89.262080389218,
            179.99220798277538,
            20003925.854,
        ),
        (
            89.333123580033,
            -89.333123580033,
            179.99295812360148,
            20003926.881,
        ),
        (
            56.320923501171,
            -56.320923501171,
            179.664_747_671_772_9,
            19993558.287,
        ),
        (
            52.784459512564,
            -52.78445951256399,
            179.634_407_464_943_8,
            19991596.095,
        ),
        (
            48.522876735459,
            -48.52287673545898,
            179.59972045622308,
            19989144.774,
        ),
        (45.0, -45.0, 179.572719, 19987083.007),
    ] {
        close(
            ground_distance_m(point(lat1, 0.0), point(lat2, lon2)),
            distance,
            0.0005,
        );
    }
    // GeodSolve78: published failure-to-converge case for the NGS calculator.
    let inverse = SurfaceGeodesic::wgs84().inverse(point(27.2, 0.0), point(-27.1, 179.5));
    close(inverse.distance_m, 19974354.765767, 0.5e-6);
    close(inverse.initial_azimuth_deg, 45.82468716758, 0.5e-11);
    close(inverse.final_azimuth_deg, 134.22776532670, 0.5e-11);
    // GeodSolve33: canonical equatorial antipode and a near-antipode.
    for (lon, distance, heading, end_heading) in [
        (179.5, 19980862.0, 55.96650, 124.03350),
        (180.0, 20003931.0, 0.0, 180.0),
    ] {
        let inverse = SurfaceGeodesic::wgs84().inverse(point(0.0, 0.0), point(0.0, lon));
        close(inverse.distance_m, distance, 0.5);
        close(
            angle_error(inverse.initial_azimuth_deg, heading),
            0.0,
            0.5e-5,
        );
        close(
            angle_error(inverse.final_azimuth_deg, end_heading),
            0.0,
            0.5e-5,
        );
    }
    // GeodSolve73: backwards from the north pole.
    let direct = SurfaceGeodesic::wgs84()
        .direct(point(90.0, 10.0), 180.0, -1e6)
        .unwrap();
    close(direct.point.latitude_deg(), 81.04623, 0.5e-5);
    close(
        angle_error(direct.point.longitude_deg(), -170.0),
        0.0,
        0.5e-5,
    );
    close(direct.final_azimuth_deg, 0.0, 0.5e-5);
}

#[test]
fn shortest_paths_and_midpoints_across_global_routes() {
    let geodesic = SurfaceGeodesic::wgs84();
    for (start, end) in [
        (point(51.5, -0.12), point(-41.3, 174.78)), // Europe to New Zealand
        (point(52.5, 13.4), point(-33.9, 151.2)),   // Europe to Australia
        (point(10.0, 179.8), point(10.1, -179.7)),
        (point(0.0, 0.0), point(0.0, 180.0)),
        (point(89.0, 0.0), point(89.0, 180.0)),
        (point(-90.0, 13.0), point(90.0, 120.0)),
        (point(51.5, -0.12), point(51.5, -0.12)),
    ] {
        let path = geodesic.path(start.with_height(200.0).unwrap(), end);
        let total = path.inverse().distance_m;
        let midpoint = path.midpoint();
        assert_eq!(midpoint.height_m(), 0.0);
        close(ground_distance_m(start, midpoint), total / 2.0, 1e-7);
        close(ground_distance_m(midpoint, end), total / 2.0, 1e-7);
        close(
            ground_distance_m(path.point_at_fraction(1.0).unwrap(), end),
            0.0,
            1e-7,
        );
        close(
            ground_distance_m(start, end),
            ground_distance_m(end, start),
            1e-7,
        );
        for fraction in [0.0, 0.1, 0.75, 1.0] {
            let along = path.point_at_fraction(fraction).unwrap();
            close(ground_distance_m(start, along), total * fraction, 1e-7);
        }
    }
    // A high-latitude surface midpoint bows poleward and stays at surface height.
    let start = point(70.0, -60.0);
    let end = point(70.0, 60.0);
    let surface = geodesic.path(start, end).midpoint();
    assert!(surface.latitude_deg() > 79.0);
    let chord = ((start.to_ecef() + end.to_ecef()) / 2.0)
        .to_geodetic()
        .unwrap();
    assert!(chord.height_m() < -100_000.0);
}

#[test]
fn spherical_compatibility_and_planetary_geodesics() {
    close(
        great_circle_m(point(0.0, 0.0), point(0.0, 90.0)),
        MEAN_EARTH_RADIUS_M * std::f64::consts::FRAC_PI_2,
        1e-8,
    );
    close(
        great_circle_m(point(0.0, 0.0), point(0.0, 180.0)),
        MEAN_EARTH_RADIUS_M * std::f64::consts::PI,
        1e-8,
    );
    close(
        great_circle_m(point(90.0, 0.0), point(90.0, 180.0)),
        0.0,
        1e-8,
    );
    close(
        great_circle_m(point(0.0, 0.0), point(0.0, 1e-12)),
        MEAN_EARTH_RADIUS_M * 1e-12_f64.to_radians(),
        1e-18,
    );
    close(
        ground_distance_m(point(0.0, 0.0), point(0.0, 1e-12)),
        WGS84_A * 1e-12_f64.to_radians(),
        1e-11,
    );
    let lunar_radius = 1_737_400.0;
    let sphere = SurfaceGeodesic::new(Ellipsoid::new(lunar_radius, 0.0).unwrap()).unwrap();
    close(
        sphere.inverse(point(0.0, 0.0), point(0.0, 90.0)).distance_m,
        lunar_radius * std::f64::consts::FRAC_PI_2,
        1e-8,
    );
    let mars = SurfaceGeodesic::new(Ellipsoid::new(3_396_190.0, 1.0 / 169.8).unwrap()).unwrap();
    let inverse = mars.inverse(point(45.0, 70.0), point(-44.9, -109.9));
    let direct = mars
        .direct(
            point(45.0, 70.0),
            inverse.initial_azimuth_deg,
            inverse.distance_m,
        )
        .unwrap();
    close(
        mars.inverse(direct.point, point(-44.9, -109.9)).distance_m,
        0.0,
        1e-7,
    );
}

#[test]
fn geodesic_domain_checks() {
    for ellipsoid in [
        Ellipsoid::new(6e6, 0.1).unwrap(),
        Ellipsoid::new(1e200, 0.0).unwrap(),
        Ellipsoid::new(1e-200, 0.0).unwrap(),
    ] {
        assert!(SurfaceGeodesic::new(ellipsoid).is_err());
    }
    let solver = SurfaceGeodesic::wgs84();
    let start = point(0.0, 0.0);
    for (azimuth, distance) in [(f64::NAN, 0.0), (0.0, f64::INFINITY)] {
        assert!(solver.direct(start, azimuth, distance).is_err());
    }
    let path = solver.path(start, point(0.0, 1.0));
    for fraction in [-0.1, 1.1, f64::NAN] {
        assert!(path.point_at_fraction(fraction).is_err());
    }
    assert!(path.position(-1.0).is_err());
    assert!(path.position(path.inverse().distance_m + 1.0).is_err());
}

#[test]
fn curvature_equator_pole_sphere_and_euler() {
    let equator = point(0.0, 0.0).curvature_radii();
    close(equator.prime_vertical_m(), WGS84_A, 1e-9);
    close(equator.meridional_m(), 6_335_439.327292819, 1e-8);
    let pole = point(90.0, 123.0).curvature_radii();
    close(pole.meridional_m(), WGS84_A / (1.0 - WGS84_F), 1e-8);
    close(pole.prime_vertical_m(), pole.meridional_m(), 1e-8);
    let latitude = point(52.0, 5.0);
    let radii = latitude.curvature_radii();
    close(
        radii.along_azimuth_m(0.0).unwrap(),
        radii.meridional_m(),
        1e-8,
    );
    close(
        radii.along_azimuth_m(90.0).unwrap(),
        radii.prime_vertical_m(),
        1e-8,
    );
    close(
        radii.along_azimuth_m(45.0).unwrap(),
        2.0 / (1.0 / radii.meridional_m() + 1.0 / radii.prime_vertical_m()),
        1e-8,
    );
    close(
        Ellipsoid::WGS84
            .radius_along_azimuth_m(latitude, 225.0)
            .unwrap(),
        radii.along_azimuth_m(45.0).unwrap(),
        1e-8,
    );
    assert!(radii.along_azimuth_m(f64::NAN).is_err());
    let sphere = Ellipsoid::new(1_737_400.0, 0.0).unwrap();
    for lat in [-90.0, -45.0, 0.0, 45.0, 90.0] {
        let radii = sphere.curvature_radii(point(lat, 0.0)).unwrap();
        close(radii.meridional_m(), 1_737_400.0, 1e-8);
        close(radii.prime_vertical_m(), 1_737_400.0, 1e-8);
        close(radii.along_azimuth_m(37.0).unwrap(), 1_737_400.0, 1e-8);
    }
    // An independent differential ECEF check: M and N cos(lat) scale derivatives.
    let delta = 1e-5;
    let north = point(52.0 + delta, 5.0)
        .to_ecef()
        .distance(point(52.0 - delta, 5.0).to_ecef());
    let east = point(52.0, 5.0 + delta)
        .to_ecef()
        .distance(point(52.0, 5.0 - delta).to_ecef());
    close(
        north / (2.0 * delta.to_radians()),
        radii.meridional_m(),
        0.01,
    );
    close(
        east / (2.0 * delta.to_radians()),
        radii.prime_vertical_m() * 52.0_f64.to_radians().cos(),
        0.01,
    );
}

#[test]
fn cartesian_interop_and_validated_construction() {
    let original = [1.0, 2.0, 3.0];
    let ecef: Ecef = original.into();
    let vector: Vec3 = ecef.into();
    assert_eq!(<[f64; 3]>::from(vector), original);
    assert_eq!(<[f64; 3]>::from(Ecef::from(Vec3::from(original))), original);
    let displacement = Vec3::from([3.0, 4.0, 5.0]);
    assert_eq!((ecef + displacement) - ecef, displacement);
    assert_eq!((ecef + displacement) - displacement, ecef);
    assert_eq!((ecef + ecef) / 2.0, ecef);
    assert_eq!(ecef * 2.0, ecef + ecef);
    assert_eq!(
        Geodetic::from_validated(10.0, 370.0, -25.0),
        Geodetic::new(10.0, 10.0, -25.0).unwrap()
    );
    assert!(std::panic::catch_unwind(|| Geodetic::from_validated(91.0, 0.0, 0.0)).is_err());
}
