use geoid_atlas::Error;
use geoid_atlas::coordinates::*;
use geoid_atlas::frames::*;
use geoid_atlas::gravity::*;
use geoid_atlas::interferometry::*;
use geoid_atlas::raster::*;
use geoid_atlas::rf::*;

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a} differs from {b}, tolerance {tol}"
    );
}
fn vector(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn epoch(t: f64) -> Epoch {
    Epoch::new(t).unwrap()
}

#[test]
fn ecef_control_points_and_round_trips() {
    let equator = Geodetic::new(0.0, 0.0, 0.0).unwrap().to_ecef();
    close(equator.x, WGS84_A, 1e-9);
    close(equator.y, 0.0, 1e-9);
    close(equator.z, 0.0, 1e-9);
    for lat in [-90.0, -89.999, -45.0, 0.0, 40.0, 89.999, 90.0] {
        for lon in [-179.99, -105.0, 0.0, 179.99] {
            for h in [-2000.0, 0.0, 400_000.0, 1e10] {
                let p = Geodetic::new(lat, lon, h).unwrap();
                let q = p.to_ecef().to_geodetic().unwrap();
                close(p.latitude_deg(), q.latitude_deg(), 1e-9);
                if lat.abs() < 90.0 {
                    close(p.longitude_deg(), q.longitude_deg(), 1e-9);
                }
                close(p.height_m(), q.height_m(), 3e-6);
            }
        }
    }
    assert!(
        Ecef {
            x: 0.0,
            y: 0.0,
            z: 0.0
        }
        .to_geodetic()
        .is_err()
    );
    assert!(Geodetic::new(f64::NAN, 0.0, 0.0).is_err());
}

#[test]
fn planetary_ellipsoid_and_subsurface_round_trip() {
    let mars = Ellipsoid::new(3_396_190.0, 1.0 / 169.8).unwrap();
    let p = Geodetic::new(-35.0, 120.0, -5000.0).unwrap();
    let q = mars.to_geodetic(mars.to_cartesian(p)).unwrap();
    close(q.latitude_deg(), p.latitude_deg(), 1e-10);
    close(q.longitude_deg(), p.longitude_deg(), 1e-10);
    close(q.height_m(), -5000.0, 1e-6);
    let sphere = Ellipsoid::new(1000.0, 0.0).unwrap();
    close(
        sphere
            .to_cartesian(Geodetic::new(90.0, 0.0, 0.0).unwrap())
            .z,
        1000.0,
        1e-10,
    );
    assert!(Ellipsoid::new(-1.0, 0.0).is_err());
}

#[test]
fn utm_external_control_and_southern_round_trip() {
    // EPSG:32613 control point, WGS84 40N 105W (central meridian).
    let zone = Utm::new(13, Hemisphere::North).unwrap();
    let p = Geodetic::new(40.0, -105.0, 123.0).unwrap();
    let (e, n) = zone.project(p).unwrap();
    close(e, 500_000.0, 1e-6);
    close(n, 4_427_757.218_738, 0.001);
    let q = zone.unproject(e, n, 123.0).unwrap();
    close(q.latitude_deg(), 40.0, 1e-7);
    close(q.longitude_deg(), -105.0, 1e-7);
    let south = Utm::new(56, Hemisphere::South).unwrap();
    let p = Geodetic::new(-33.86, 151.21, 50.0).unwrap();
    let (e, n) = south.project(p).unwrap();
    let q = south.unproject(e, n, 50.0).unwrap();
    close(q.latitude_deg(), p.latitude_deg(), 1e-7);
    close(q.longitude_deg(), p.longitude_deg(), 1e-7);
    assert!(zone.project(p).is_err());
    assert!(Utm::new(0, Hemisphere::North).is_err());
}

#[test]
fn mercator_and_local_enu() {
    let origin = Geodetic::new(40.0, -105.0, 1600.0).unwrap();
    let (x, y) = Crs::WebMercator.project(origin).unwrap();
    let q = Crs::WebMercator.unproject(x, y, 1600.0).unwrap();
    close(q.latitude_deg(), 40.0, 1e-10);
    close(q.longitude_deg(), -105.0, 1e-10);
    assert!(
        Crs::WebMercator
            .project(Geodetic::new(90.0, 0.0, 0.0).unwrap())
            .is_err()
    );
    let frame = LocalFrame::new(origin);
    let enu = Enu {
        east: 100.0,
        north: -50.0,
        up: 10.0,
    };
    let round = frame.to_enu(frame.to_ecef(enu));
    close(round.east, enu.east, 1e-8);
    close(round.north, enu.north, 1e-8);
    close(round.up, enu.up, 1e-8);
}

#[test]
fn rotating_frame_velocity_and_inverse() {
    let mut graph = FrameGraph::new("inertial");
    let root = graph.root();
    let rotating = graph
        .add(
            "planet",
            root,
            LinearMotion {
                reference_epoch: epoch(0.0),
                origin: State {
                    position_m: vector(1e8, 2.0, 3.0),
                    velocity_m_s: vector(10.0, 0.0, 0.0),
                },
                initial_rotation: Rotation::IDENTITY,
                angular_velocity_rad_s: vector(0.0, 0.0, 0.01),
            },
        )
        .unwrap();
    let local = State::stationary(vector(100.0, 0.0, 0.0));
    let state = graph.transform(local, rotating, root, epoch(0.0)).unwrap();
    close(state.velocity_m_s.x, 10.0, 1e-12);
    close(state.velocity_m_s.y, 1.0, 1e-12);
    let back = graph.transform(state, root, rotating, epoch(0.0)).unwrap();
    close((back.position_m - local.position_m).norm(), 0.0, 1e-8);
    close(back.velocity_m_s.norm(), 0.0, 1e-12);
    let before = graph
        .transform(local, rotating, root, epoch(-0.01))
        .unwrap();
    let after = graph.transform(local, rotating, root, epoch(0.01)).unwrap();
    close(
        ((after.position_m - before.position_m) * (1.0 / 0.02) - state.velocity_m_s).norm(),
        0.0,
        2e-6,
    );
}

#[test]
fn local_precision_at_interstellar_translation_and_graph_ownership() {
    let mut graph = FrameGraph::new("root");
    let root = graph.root();
    let star = graph
        .add(
            "star",
            root,
            Pose {
                origin: State::stationary(vector(1e20, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let a = graph
        .add(
            "a",
            star,
            Pose {
                origin: State::stationary(vector(100.0, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let b = graph
        .add(
            "b",
            star,
            Pose {
                origin: State::stationary(vector(101.0, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let result = graph
        .transform(State::stationary(vector(0.001, 0.0, 0.0)), a, b, epoch(0.0))
        .unwrap();
    close(result.position_m.x, -0.999, 1e-12);
    let other = FrameGraph::new("other");
    assert!(
        graph
            .transform(
                State::stationary(Vec3::ZERO),
                other.root(),
                root,
                epoch(0.0)
            )
            .is_err()
    );
    assert!(
        graph
            .transform(
                State::stationary(Vec3::ZERO),
                root,
                other.root(),
                epoch(0.0)
            )
            .is_err()
    );
    assert!(Rotation::new([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).is_err());
}

#[test]
fn kepler_orbit_conserves_energy_and_returns_after_period() {
    let a: f64 = 7e6;
    let mu: f64 = 3.986_004_418e14;
    for e in [0.0, 0.4, 0.99] {
        let orbit = KeplerOrbit {
            reference_epoch: epoch(0.0),
            semi_major_axis_m: a,
            eccentricity: e,
            gravitational_parameter_m3_s2: mu,
            mean_anomaly_rad: 0.2,
            orbital_plane: Rotation::IDENTITY,
        };
        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / mu).sqrt();
        let start = orbit.pose(epoch(0.0)).unwrap().origin;
        let end = orbit.pose(epoch(period)).unwrap().origin;
        close((start.position_m - end.position_m).norm(), 0.0, 1e-6);
        for t in [0.0, 100.0, 2000.0] {
            let s = orbit.pose(epoch(t)).unwrap().origin;
            close(
                0.5 * s.velocity_m_s.dot(s.velocity_m_s) - mu / s.position_m.norm(),
                -mu / (2.0 * a),
                1e-5,
            );
        }
    }
}

#[test]
fn bilinear_rotated_raster_and_nodata() {
    let affine = Affine::new(10.0, 20.0, 2.0, 1.0, -1.0, 2.0).unwrap();
    let (x, y) = affine.world(0.25, 0.75);
    let raster = Raster::new(
        2,
        2,
        affine,
        Crs::Planetographic,
        vec![Some(0.0), Some(10.0), Some(20.0), Some(30.0)],
    )
    .unwrap();
    close(raster.bilinear_xy(x, y).unwrap(), 17.5, 1e-12);
    let missing = Raster::new(
        2,
        2,
        affine,
        Crs::Planetographic,
        vec![Some(0.0), None, Some(20.0), Some(30.0)],
    )
    .unwrap();
    assert_eq!(missing.bilinear_xy(x, y), Err(Error::NoData));
    close(missing.bilinear_xy(10.0, 20.0).unwrap(), 0.0, 1e-12);
    assert_eq!(
        raster.bilinear_xy(-100.0, 20.0),
        Err(Error::OutsideCoverage)
    );
    assert!(Affine::new(0.0, 0.0, 1.0, 1.0, 2.0, 2.0).is_err());
}

#[test]
fn ascii_dem_row_order_corners_and_missing_data() {
    let text = "ncols 2\nnrows 2\nxllcorner 10\nyllcorner 20\ncellsize 1\nNODATA_value -9999\n100 200\n300 -9999\n";
    let dem = Dem::from_esri_ascii(text, Crs::Planetographic, VerticalDatum::Orthometric).unwrap();
    close(dem.raster.bilinear_xy(10.5, 21.5).unwrap(), 100.0, 1e-12);
    close(dem.raster.bilinear_xy(10.5, 20.5).unwrap(), 300.0, 1e-12);
    assert_eq!(dem.raster.bilinear_xy(11.5, 20.5), Err(Error::NoData));
    assert!(
        Dem::from_esri_ascii(
            &text.replace("ncols 2", "ncols 3"),
            Crs::Planetographic,
            VerticalDatum::Ellipsoid
        )
        .is_err()
    );
    assert!(
        Dem::from_esri_ascii(
            &text.replace("cellsize 1", "cellsize -1"),
            Crs::Planetographic,
            VerticalDatum::Ellipsoid
        )
        .is_err()
    );
}

#[test]
fn gravity_controls_and_mass_superposition() {
    close(wgs84_normal_gravity(0.0).unwrap(), 9.780_325_335_9, 1e-12);
    close(wgs84_normal_gravity(90.0).unwrap(), 9.832_184_937_8, 1e-10);
    close(ellipsoidal_height(100.0, -30.0).unwrap(), 70.0, 1e-12);
    close(orthometric_height(70.0, -30.0).unwrap(), 100.0, 1e-12);
    let field = point_mass_field(
        Vec3::ZERO,
        &[
            PointMass {
                position_m: vector(10.0, 0.0, 0.0),
                gravitational_parameter_m3_s2: 100.0,
            },
            PointMass {
                position_m: vector(-10.0, 0.0, 0.0),
                gravitational_parameter_m3_s2: 100.0,
            },
        ],
    )
    .unwrap();
    close(field.acceleration_m_s2.norm(), 0.0, 1e-12);
    close(field.potential_m2_s2, -20.0, 1e-12);
}

#[test]
fn geometric_delay_sign_and_far_field_limit() {
    let receiver1 = State::stationary(Vec3::ZERO);
    let receiver2 = State::stationary(vector(300.0, 100.0, 0.0));
    let far = far_field(receiver1, receiver2, vector(1.0, 0.0, 0.0)).unwrap();
    close(far.seconds, -300.0 / SPEED_OF_LIGHT_M_S, 1e-20);
    let near = near_field_static(
        vector(1e15, 0.0, 0.0),
        receiver1.position_m,
        receiver2.position_m,
    )
    .unwrap();
    close(near, far.seconds, 1e-17);
    close(far.uvw_m.norm(), receiver2.position_m.norm(), 1e-10);
    let midpoint = vector(150.0, 50.0, 0.0);
    close(
        near_field_static(midpoint, receiver1.position_m, receiver2.position_m).unwrap(),
        0.0,
        1e-20,
    );
}

#[test]
fn moving_receiver_light_time_matches_analytic_solution() {
    let graph = FrameGraph::new("inertial");
    let speed = 1000.0;
    let initial_distance = 3e8;
    let receiver = Receiver {
        frame: graph.root(),
        reference_epoch: epoch(0.0),
        local_state: State {
            position_m: vector(initial_distance, 0.0, 0.0),
            velocity_m_s: vector(speed, 0.0, 0.0),
        },
    };
    let solved = arrival(&graph, Vec3::ZERO, epoch(0.0), receiver, 1e-13).unwrap();
    close(
        solved.travel_time_s,
        initial_distance / (SPEED_OF_LIGHT_M_S - speed),
        1e-12,
    );
    assert!(solved.iterations > 1);
    assert!(
        arrival(
            &graph,
            Vec3::ZERO,
            epoch(0.0),
            Receiver {
                local_state: State {
                    velocity_m_s: vector(SPEED_OF_LIGHT_M_S, 0.0, 0.0),
                    ..receiver.local_state
                },
                ..receiver
            },
            1e-9
        )
        .is_err()
    );
}

#[test]
fn medium_clocks_phase_and_expected_coherence() {
    close(
        medium_delay(100.0, 2.0).unwrap(),
        100.0 / SPEED_OF_LIGHT_M_S,
        1e-20,
    );
    let clock = ClockModel {
        reference_epoch: epoch(10.0),
        offset_s: 1e-6,
        drift_s_per_s: 1e-9,
        drift_rate_s_per_s2: 0.0,
    };
    close(clock.offset_at(epoch(20.0)).unwrap(), 1.01e-6, 1e-20);
    let phase = DelayBudget {
        geometric_s: 0.25e-9,
        ..Default::default()
    }
    .phase_rad(1e9)
    .unwrap();
    close(phase, std::f64::consts::FRAC_PI_2, 1e-12);
    let coherence = expected_coherence(0.25e-9, 0.0, 1e9).unwrap();
    close(coherence.real, 0.0, 1e-12);
    close(coherence.imaginary, 1.0, 1e-12);
    let smeared = expected_coherence(0.0, 1e-9, 1e9).unwrap();
    assert!(smeared.real < 1e-8);
}

#[test]
fn gravitational_delay_and_occultation() {
    let a = vector(1e8, 1e8, 0.0);
    let b = vector(2e8, 1e8, 0.0);
    let delay = shapiro_delay(a, b, Vec3::ZERO, 3.986e14).unwrap();
    assert!(delay > 0.0);
    close(
        delay,
        shapiro_delay(b, a, Vec3::ZERO, 3.986e14).unwrap(),
        1e-20,
    );
    assert!(
        shapiro_delay(
            vector(-1.0, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
            Vec3::ZERO,
            1.0
        )
        .is_err()
    );
    assert!(
        ray_intersects_sphere(
            vector(-2.0, 0.0, 0.0),
            vector(2.0, 0.0, 0.0),
            Vec3::ZERO,
            1.0
        )
        .unwrap()
    );
    assert!(
        !ray_intersects_sphere(
            vector(-2.0, 2.0, 0.0),
            vector(2.0, 2.0, 0.0),
            Vec3::ZERO,
            1.0
        )
        .unwrap()
    );
}

#[test]
fn relativistic_doppler_control() {
    let beta = 0.1;
    close(
        doppler_ratio(
            Vec3::ZERO,
            vector(beta * SPEED_OF_LIGHT_M_S, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
        )
        .unwrap(),
        ((1.0 - beta) / (1.0 + beta)).sqrt(),
        1e-12,
    );
    close(
        doppler_ratio(
            vector(1000.0, 0.0, 0.0),
            vector(1000.0, 0.0, 0.0),
            vector(1.0, 0.0, 0.0),
        )
        .unwrap(),
        1.0,
        1e-12,
    );
}

#[test]
fn rf_friis_lognormal_and_independent_power() {
    close(
        free_space_path_loss_db(1000.0, 1e9).unwrap(),
        92.447_783_221_9,
        1e-9,
    );
    let power = LinkBudget {
        transmitted_power_dbm: 30.0,
        transmit_gain_dbi: 0.0,
        receive_gain_dbi: 0.0,
        additional_loss_db: 0.0,
        shadowing_std_db: 0.0,
    }
    .predict(1000.0, 1e9)
    .unwrap();
    close(power.mean_dbm(), 30.0 - 92.447_783_221_9, 1e-9);
    let uncertain = PowerDistribution::new(0.0, 3.0).unwrap();
    assert!(uncertain.moments().unwrap().mean_w > uncertain.median_w());
    let total = independent_power_sum(&[uncertain, uncertain]).unwrap();
    close(
        total.mean_w,
        2.0 * uncertain.moments().unwrap().mean_w,
        1e-14,
    );
    close(
        total.variance_w2,
        2.0 * uncertain.moments().unwrap().variance_w2,
        1e-14,
    );
    assert!(PowerDistribution::new(0.0, -1.0).is_err());
    assert!(free_space_path_loss_db(0.0, 1e9).is_err());
}

#[test]
fn subsurface_material_vacuum_and_conductor() {
    let vacuum = Material {
        relative_permittivity: 1.0,
        relative_permeability: 1.0,
        conductivity_s_m: 0.0,
    }
    .propagation_constant(1e9)
    .unwrap();
    close(vacuum.attenuation_nepers_m, 0.0, 1e-14);
    close(
        vacuum.phase_rad_m,
        2.0 * std::f64::consts::PI * 1e9 / SPEED_OF_LIGHT_M_S,
        1e-12,
    );
    let rock = Material {
        relative_permittivity: 9.0,
        relative_permeability: 1.0,
        conductivity_s_m: 0.01,
    }
    .propagation_constant(1e8)
    .unwrap();
    assert!(rock.attenuation_nepers_m > 0.0);
    assert!(rock.absorption_loss_db(100.0).unwrap() > 0.0);
    let metal = Material {
        relative_permittivity: 1.0,
        relative_permeability: 1.0,
        conductivity_s_m: 1e7,
    }
    .propagation_constant(1e3)
    .unwrap();
    close(
        metal.attenuation_nepers_m
            / (std::f64::consts::PI * 1e3 * 1.256_637_061_27e-6 * 1e7).sqrt(),
        1.0,
        1e-12,
    );
}

#[test]
fn weather_and_plasma_corrections_keep_group_and_phase_signs() {
    let n = tropospheric_refractive_index(293.0, 101_325.0, 1000.0).unwrap();
    assert!(n > 1.0002 && n < 1.0004);
    assert!(tropospheric_refractive_index(293.0, 10.0, 100.0).is_err());
    let plasma = ionospheric_delay(1e17, 1e9).unwrap();
    close(plasma.group_delay_s, 4.03 / SPEED_OF_LIGHT_M_S, 1e-20);
    close(plasma.phase_delay_s, -plasma.group_delay_s, 1e-20);
}
