use geoid_atlas::coordinates::{Ellipsoid, Geodetic};
use geoid_atlas::estimation::Covariance6;
use geoid_atlas::formats::{Helmert, read_world_file};
use geoid_atlas::frames::*;
use geoid_atlas::gravity::GeoidModel;
use geoid_atlas::reference::*;
use geoid_atlas::units::*;
use geoid_atlas::{Error, Result};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn unit_boundaries_distinguish_feet_angles_and_weather_units() {
    close(
        Length::new(1.0, LengthUnit::AstronomicalUnit)
            .unwrap()
            .metres(),
        149_597_870_700.0,
        0.0,
    );
    let survey = Length::new(1_000_000.0, LengthUnit::UsSurveyFoot)
        .unwrap()
        .metres();
    let international = Length::new(1_000_000.0, LengthUnit::InternationalFoot)
        .unwrap()
        .metres();
    close(survey - international, 0.6096012192, 1e-10);
    close(
        Angle::new(6.0, AngleUnit::HourAngle)
            .unwrap()
            .in_unit(AngleUnit::Degree)
            .unwrap(),
        90.0,
        1e-12,
    );
    close(
        Duration::new(1.0, DurationUnit::JulianYear)
            .unwrap()
            .seconds(),
        31_557_600.0,
        0.0,
    );
    close(
        Frequency::new(1.57542, FrequencyUnit::Gigahertz)
            .unwrap()
            .hertz(),
        1_575_420_000.0,
        1e-6,
    );
    close(
        Pressure::new(1013.25, PressureUnit::Hectopascal)
            .unwrap()
            .pascals(),
        101_325.0,
        0.0,
    );
    close(
        MagneticFluxDensity::new(50_000.0, MagneticUnit::Nanotesla)
            .unwrap()
            .teslas(),
        50e-6,
        1e-20,
    );
    close(
        Temperature::new(32.0, TemperatureUnit::Fahrenheit)
            .unwrap()
            .kelvin(),
        273.15,
        0.0,
    );
    assert!(Temperature::new(-1.0, TemperatureUnit::Kelvin).is_err());
    assert!(
        Temperature::new(f64::MAX, TemperatureUnit::Kelvin)
            .unwrap()
            .in_unit(TemperatureUnit::Fahrenheit)
            .is_err()
    );
    assert!(Length::new(f64::MAX, LengthUnit::Parsec).is_err());
    assert!(Angle::new(f64::NAN, AngleUnit::Degree).is_err());
}

struct Geoid(f64);
impl GeoidModel for Geoid {
    fn undulation_m(&self, _: Geodetic) -> Result<f64> {
        Ok(self.0)
    }
}
struct Missing;
impl GeoidModel for Missing {
    fn undulation_m(&self, _: Geodetic) -> Result<f64> {
        Err(Error::NoData)
    }
}

#[test]
fn vertical_reference_change_and_horizontal_shift_are_separate() {
    let mut graph = FrameGraph::new("datum A");
    let root = graph.root();
    let shifted = graph
        .add(
            "datum B",
            root,
            Pose {
                origin: State::stationary(v(3.0, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let ga = Geoid(30.0);
    let gb = Geoid(20.0);
    let a = GeodeticReference {
        frame: root,
        ellipsoid: Ellipsoid::WGS84,
        height: HeightReference::Orthometric(&ga),
    };
    let b = GeodeticReference {
        frame: shifted,
        ellipsoid: Ellipsoid::WGS84,
        height: HeightReference::Orthometric(&gb),
    };
    let point = GeographicCoordinate {
        latitude_deg: 0.0,
        longitude_deg: 0.0,
        height_m: 100.0,
    };
    let converted = a
        .convert(point, &b, &graph, Epoch::new(0.0).unwrap())
        .unwrap();
    close(converted.height_m, 107.0, 1e-8); // hA=130; xB=xA-3; HB=127-20
    close(
        b.convert(converted, &a, &graph, Epoch::new(0.0).unwrap())
            .unwrap()
            .height_m,
        100.0,
        1e-8,
    );
    let missing = GeodeticReference {
        frame: root,
        ellipsoid: Ellipsoid::WGS84,
        height: HeightReference::Orthometric(&Missing),
    };
    assert_eq!(missing.to_cartesian(point), Err(Error::NoData));
    let other = FrameGraph::new("unrelated");
    assert!(
        a.convert(point, &b, &other, Epoch::new(0.0).unwrap())
            .is_err()
    );
}

#[test]
fn lunar_surface_event_is_not_relocated_by_earth_conversion() {
    let mut graph = FrameGraph::new("Earth fixed illustrative instant");
    let root = graph.root();
    let lunar = graph
        .add(
            "Moon fixed",
            root,
            Pose {
                origin: State::stationary(v(384_400_000.0, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let moon = GeodeticReference {
        frame: lunar,
        ellipsoid: Ellipsoid::new(1_737_400.0, 0.0).unwrap(),
        height: HeightReference::Ellipsoidal,
    };
    let earth = GeodeticReference {
        frame: root,
        ellipsoid: Ellipsoid::WGS84,
        height: HeightReference::Ellipsoidal,
    };
    let surface = GeographicCoordinate {
        latitude_deg: 0.0,
        longitude_deg: 0.0,
        height_m: 0.0,
    };
    let converted = moon
        .convert(surface, &earth, &graph, Epoch::new(0.0).unwrap())
        .unwrap();
    close(
        converted.height_m,
        384_400_000.0 + 1_737_400.0 - 6_378_137.0,
        1e-6,
    );
    let round = earth
        .convert(converted, &moon, &graph, Epoch::new(0.0).unwrap())
        .unwrap();
    close(round.height_m, 0.0, 1e-6);
}

#[test]
fn mars_tangent_frame_and_centric_latitude() {
    let mars = Ellipsoid::new(3_396_190.0, 1.0 / 169.8).unwrap();
    let site = Geodetic::new(45.0, 20.0, -2000.0).unwrap();
    let pose = tangent_pose(mars, site).unwrap();
    let normal = pose.rotation.apply(v(0.0, 0.0, 1.0));
    close(normal.z, 0.5_f64.sqrt(), 1e-14);
    let centric = Planetocentric::from_cartesian(pose.origin.position_m).unwrap();
    assert!(centric.latitude_deg < 45.0);
    assert!((centric.to_cartesian().unwrap() - pose.origin.position_m).norm() < 1e-8);
    assert!(Planetocentric::from_cartesian(Vec3::ZERO).is_err());
    for latitude in [-90.0, 90.0] {
        assert!(tangent_pose(mars, Geodetic::new(latitude, 0.0, 0.0).unwrap()).is_ok());
    }
}

#[test]
fn helmert_wgs72_to_wgs84_published_rounded_control() {
    // EPSG position-vector example, rounded output coordinates (centimetre precision).
    let h =
        Helmert::from_proj("+proj=helmert +z=4.5 +rz=0.554 +s=0.219 +convention=position_vector")
            .unwrap();
    let result = h
        .at_decimal_year(2000.0)
        .unwrap()
        .apply(State::stationary(v(3_657_660.66, 255_768.55, 5_201_382.11)))
        .unwrap();
    close(result.position_m.x, 3_657_660.78, 0.01);
    close(result.position_m.y, 255_778.43, 0.01);
    close(result.position_m.z, 5_201_387.75, 0.01);
}

#[test]
fn dynamic_helmert_velocity_covariance_and_exact_inverse() {
    let h = Helmert::from_proj("+proj=helmert +x=1 +y=-2 +z=3 +rx=0.1 +ry=-0.2 +rz=0.3 +s=2 +dx=0.02 +dy=-0.03 +dz=0.01 +drx=0.002 +dry=-0.003 +drz=0.004 +ds=0.01 +t_epoch=2010 +convention=coordinate_frame").unwrap();
    let state = State {
        position_m: v(6e6, 2e6, -3e6),
        velocity_m_s: v(0.1, -0.2, 0.3),
    };
    let transform = h.at_decimal_year(2020.0).unwrap();
    let result = transform.apply(state).unwrap();
    let round = transform.unapply(result).unwrap();
    assert!((round.position_m - state.position_m).norm() < 2e-9);
    assert!((round.velocity_m_s - state.velocity_m_s).norm() < 2e-16);
    let dt = 100.0;
    let before = h
        .at_decimal_year(2020.0 - dt / 31_557_600.0)
        .unwrap()
        .apply(State {
            position_m: state.position_m - state.velocity_m_s * dt,
            ..state
        })
        .unwrap();
    let after = h
        .at_decimal_year(2020.0 + dt / 31_557_600.0)
        .unwrap()
        .apply(State {
            position_m: state.position_m + state.velocity_m_s * dt,
            ..state
        })
        .unwrap();
    assert!(
        ((after.position_m - before.position_m) * (0.5 / dt) - result.velocity_m_s).norm() < 1e-11
    );
    let covariance = transform
        .covariance(Covariance6::isotropic(1.0, 0.01).unwrap())
        .unwrap();
    assert!(covariance.matrix()[0][3].abs() > 0.0); // scale/rotation-rate coupling preserved
}

#[test]
fn helmert_conventions_and_strict_format_errors() {
    let position = Helmert::from_proj("+proj=helmert +rz=1 +convention=position_vector").unwrap();
    let coordinate =
        Helmert::from_proj("+proj=helmert +rz=1 +convention=coordinate_frame").unwrap();
    let state = State::stationary(v(1e6, 0.0, 0.0));
    close(
        position
            .at_decimal_year(2000.0)
            .unwrap()
            .apply(state)
            .unwrap()
            .position_m
            .y,
        -coordinate
            .at_decimal_year(2000.0)
            .unwrap()
            .apply(state)
            .unwrap()
            .position_m
            .y,
        0.0,
    );
    for invalid in [
        "+proj=helmert +rx=1",
        "+proj=helmert +dx=1",
        "+proj=helmert +x=NaN",
        "+proj=helmert +x=0 +x=1",
        "+proj=helmert +exact",
        "+proj=pipeline",
        "+proj=helmert +grids=foo",
        "+proj=helmert +convention=bogus",
    ] {
        assert!(Helmert::from_proj(invalid).is_err(), "{invalid}");
    }
    assert!(
        Helmert::from_proj("+proj=helmert +s=-1000000")
            .unwrap()
            .at_decimal_year(2020.0)
            .is_err()
    );
    assert!(
        Helmert::from_proj("+proj=helmert +rx=1000 +convention=position_vector")
            .unwrap()
            .at_decimal_year(2020.0)
            .is_err()
    );
}

#[test]
fn rotated_image_world_file_uses_pixel_centres() {
    let affine = read_world_file("2\n0.5\n-0.25\n-3\n100\n200\n").unwrap();
    assert_eq!(affine.world(0.0, 0.0), (100.0, 200.0));
    let world = affine.world(4.0, 3.0);
    let pixel = affine.pixel(world.0, world.1).unwrap();
    close(pixel.0, 4.0, 1e-14);
    close(pixel.1, 3.0, 1e-14);
    for invalid in [
        "1\n2",
        "1\n0\n0\n0\n0\n0",
        "NaN\n0\n0\n1\n0\n0",
        "1\n0\n0\n1\n0\n0\n99",
    ] {
        assert!(read_world_file(invalid).is_err());
    }
}
