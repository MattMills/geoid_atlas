use geoid_atlas::echo::*;
use geoid_atlas::environment::*;
use geoid_atlas::frames::*;
use geoid_atlas::fusion::*;
use geoid_atlas::interferometry::SPEED_OF_LIGHT_M_S as C;
use geoid_atlas::spacetime::*;
use geoid_atlas::{Error, Result};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn epoch(t: f64) -> Epoch {
    Epoch::new(t).unwrap()
}
fn world(position: Vec3, velocity: Vec3) -> ConstantAcceleration {
    ConstantAcceleration {
        reference_epoch: epoch(0.0),
        state: State {
            position_m: position,
            velocity_m_s: velocity,
        },
        acceleration_m_s2: Vec3::ZERO,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

struct UniformB(Vec3);
impl MagneticField for UniformB {
    fn teslas_at(&self, _: Vec3, _: Epoch) -> Result<Vec3> {
        Ok(self.0)
    }
}
struct UniformNe(f64);
impl ElectronDensity for UniformNe {
    fn per_m3_at(&self, _: Vec3, _: Epoch) -> Result<f64> {
        Ok(self.0)
    }
}

#[test]
fn polar_field_axes_and_dipole_symmetry() {
    let dipole = DipoleField {
        centre_m: Vec3::ZERO,
        moment_a_m2: v(0.0, 0.0, 8e22),
        minimum_radius_m: 6e6,
    };
    let north = dipole.teslas_at(v(0.0, 0.0, 6.4e6), epoch(0.0)).unwrap();
    let south = dipole.teslas_at(v(0.0, 0.0, -6.4e6), epoch(0.0)).unwrap();
    let equator = dipole.teslas_at(v(6.4e6, 0.0, 0.0), epoch(0.0)).unwrap();
    close(north.z, south.z, 0.0);
    close(north.z, -2.0 * equator.z, 0.0);
    for b in [north, south] {
        let axes = field_aligned_basis(b, v(1.0, 0.0, 0.0)).unwrap();
        assert!((axes.apply(v(0.0, 0.0, 1.0)) - b.normalized().unwrap()).norm() < 1e-14);
    }
    assert!(field_aligned_basis(Vec3::ZERO, v(1.0, 0.0, 0.0)).is_err());
    assert!(field_aligned_basis(north, north).is_err());
    assert_eq!(
        dipole.teslas_at(Vec3::ZERO, epoch(0.0)),
        Err(Error::OutsideCoverage)
    );
}

#[test]
fn vector_field_frame_adapter_rotates_without_translating_field() {
    let mut graph = FrameGraph::new("inertial");
    let root = graph.root();
    let frame = graph
        .add(
            "sensor",
            root,
            Pose {
                origin: State::stationary(v(1000.0, 2000.0, 3000.0)),
                rotation: Rotation::axis_angle(v(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2)
                    .unwrap(),
                angular_velocity_rad_s: v(1.0, 0.0, 0.0),
            },
        )
        .unwrap();
    let field = UniformB(v(0.0, 0.0, 50e-6));
    let adapter = FramedMagneticField {
        graph: &graph,
        frame,
        field: &field,
    };
    let b = adapter
        .teslas_at(v(1000.0, 2000.0, 3000.0), epoch(0.0))
        .unwrap();
    close(b.x, 50e-6, 1e-20);
    close(b.y, 0.0, 1e-20);
    close(b.z, 0.0, 1e-20);
}

#[test]
fn rk4_field_lines_respect_analytic_dipole_invariant_and_coverage() {
    let dipole = DipoleField {
        centre_m: Vec3::ZERO,
        moment_a_m2: v(0.0, 0.0, 8e22),
        minimum_radius_m: 6e6,
    };
    // Dipole field lines satisfy r = L sin²(theta); theta is colatitude.
    let start = v(9.6e6, 0.0, 0.0);
    let line = trace_field_line(&dipole, start, epoch(0.0), 10_000.0, 100).unwrap();
    for p in line {
        close(p.norm().powi(3) / (p.x * p.x + p.y * p.y), 9.6e6, 0.01);
    }
    let uniform = UniformB(v(0.0, 0.0, 1.0));
    let north = trace_field_line(&uniform, Vec3::ZERO, epoch(0.0), 10.0, 3).unwrap();
    let south = trace_field_line(&uniform, Vec3::ZERO, epoch(0.0), -10.0, 3).unwrap();
    assert_eq!(north[3], v(0.0, 0.0, 30.0));
    assert_eq!(south[3], v(0.0, 0.0, -30.0));
    assert!(trace_field_line(&UniformB(Vec3::ZERO), start, epoch(0.0), 10.0, 1).is_err());
    assert!(trace_field_line(&dipole, v(0.0, 0.0, 6_000_001.0), epoch(0.0), -10.0, 1).is_err());
}

#[test]
fn plasma_columns_faraday_sign_and_frequency_scaling() {
    let ne = UniformNe(1e11);
    let magnetic = UniformB(v(0.0, 0.0, 50e-6));
    let start = Vec3::ZERO;
    let end = v(0.0, 0.0, 100_000.0);
    let column = integrate_plasma(&ne, &magnetic, start, end, epoch(0.0), epoch(1.0), 10).unwrap();
    close(column.electron_column_per_m2, 1e16, 1.0);
    close(column.electron_magnetic_column, 5e11, 1e-3);
    let reverse = integrate_plasma(&ne, &magnetic, end, start, epoch(0.0), epoch(1.0), 10).unwrap();
    close(
        reverse.electron_column_per_m2,
        column.electron_column_per_m2,
        0.0,
    );
    close(
        reverse.electron_magnetic_column,
        -column.electron_magnetic_column,
        0.0,
    );
    close(
        column.faraday_rotation_rad(1e8).unwrap(),
        4.0 * column.faraday_rotation_rad(2e8).unwrap(),
        1e-14,
    );
    assert!(
        integrate_plasma(
            &UniformNe(-1.0),
            &magnetic,
            start,
            end,
            epoch(0.0),
            epoch(1.0),
            10
        )
        .is_err()
    );
    assert!(integrate_plasma(&ne, &magnetic, start, end, epoch(1.0), epoch(0.0), 10).is_err());
}

#[test]
fn plasma_provider_receives_position_and_light_path_time() {
    struct Timed;
    impl ElectronDensity for Timed {
        fn per_m3_at(&self, p: Vec3, t: Epoch) -> Result<f64> {
            Ok(1e10 * (1.0 + p.z / 100_000.0 + t.duration_since(epoch(0.0))))
        }
    }
    let column = integrate_plasma(
        &Timed,
        &UniformB(Vec3::ZERO),
        Vec3::ZERO,
        v(0.0, 0.0, 100_000.0),
        epoch(0.0),
        epoch(1.0),
        2,
    )
    .unwrap();
    close(column.electron_column_per_m2, 2e15, 1.0);
}

fn zero_prediction(_: &dyn Worldline) -> Result<Vec<f64>> {
    Ok(vec![0.0, 0.0])
}

#[test]
fn correlated_gaussian_blocks_do_not_double_count_common_errors() {
    let model = GaussianObservation::new(
        vec![1.0, 1.0],
        vec![vec![1.0, 0.8], vec![0.8, 1.0]],
        zero_prediction,
    )
    .unwrap();
    let independent = GaussianObservation::new(
        vec![1.0, 1.0],
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        zero_prediction,
    )
    .unwrap();
    let target = world(Vec3::ZERO, Vec3::ZERO);
    close(model.evaluate(&target).unwrap(), 1.0 / 1.8, 1e-14);
    close(independent.evaluate(&target).unwrap(), 1.0, 1e-14);
    let mixed = GaussianObservation::new(
        vec![1.0, 1e-9],
        vec![vec![1.0, 0.8e-9], vec![0.8e-9, 1e-18]],
        zero_prediction,
    )
    .unwrap();
    close(
        mixed.evaluate(&target).unwrap(),
        model.evaluate(&target).unwrap(),
        1e-14,
    );
}

#[test]
fn gaussian_adapter_rejects_invalid_covariance_predictions_and_missing_coverage() {
    for covariance in [
        vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        vec![vec![1.0, 2.0], vec![2.0, 1.0]],
        vec![vec![1.0, 0.0], vec![0.1, 1.0]],
        vec![vec![1.0, f64::NAN], vec![0.0, 1.0]],
    ] {
        assert!(GaussianObservation::new(vec![0.0, 0.0], covariance, zero_prediction).is_err());
    }
    let target = world(Vec3::ZERO, Vec3::ZERO);
    let wrong_size = GaussianObservation::new(vec![0.0], vec![vec![1.0]], zero_prediction).unwrap();
    assert!(wrong_size.evaluate(&target).is_err());
    fn missing(_: &dyn Worldline) -> Result<Vec<f64>> {
        Err(Error::OutsideCoverage)
    }
    let absent = GaussianObservation::new(vec![0.0], vec![vec![1.0]], missing).unwrap();
    assert_eq!(absent.evaluate(&target), Err(Error::OutsideCoverage));
}

#[test]
fn backward_light_time_and_bearing_use_target_emission_event() {
    let target = world(v(1e7, 0.0, 0.0), v(1000.0, 0.0, 0.0));
    let reception = epoch(10.0);
    let emission = emission_epoch(&target, Vec3::ZERO, reception, 1e-14).unwrap();
    close(
        emission.duration_since(reception),
        -(1e7 + 1000.0 * 10.0) / (C + 1000.0),
        1e-14,
    );
    let sensor = world(Vec3::ZERO, Vec3::ZERO);
    let angle: f64 = 0.005;
    let bearing = BearingObservation {
        sensor: &sensor,
        reception_epoch: reception,
        observed_direction: v(angle.cos(), angle.sin(), 0.0),
        angular_std_rad: angle / 2.0,
        light_time_tolerance_s: 1e-14,
    };
    close(
        bearing.evaluate(&target).unwrap(),
        0.5 * (angle.tan() / (angle / 2.0)).powi(2),
        1e-12,
    );
    let behind = BearingObservation {
        observed_direction: v(-1.0, 0.0, 0.0),
        ..bearing
    };
    assert!(behind.evaluate(&target).is_err());
    let superluminal = world(v(1.0, 0.0, 0.0), v(C, 0.0, 0.0));
    assert!(emission_epoch(&superluminal, Vec3::ZERO, reception, 1e-12).is_err());
}

#[test]
fn photon_and_bearing_blocks_share_hypothesis_without_mixing_phases() {
    let sensor = world(Vec3::ZERO, Vec3::ZERO);
    let truth = world(v(1000.0, 0.0, 0.0), Vec3::ZERO);
    let wrong = world(v(1000.0, 20.0, 0.0), Vec3::ZERO);
    let bearing = BearingObservation {
        sensor: &sensor,
        reception_epoch: epoch(0.0),
        observed_direction: v(1.0, 0.0, 0.0),
        angular_std_rad: 1e-3,
        light_time_tolerance_s: 1e-12,
    };
    fn counts(w: &dyn Worldline) -> Result<Vec<f64>> {
        Ok(vec![
            100.0 / (1.0 + w.state_at(epoch(0.0))?.position_m.y.powi(2) / 400.0),
        ])
    }
    let photons = PhotonObservation {
        counts: vec![100],
        model: counts,
    };
    let mut joint = JointLikelihood::default();
    joint.add("image", "camera", bearing).unwrap();
    joint.add("xray", "counter", photons).unwrap();
    assert!(joint.evaluate(&truth).unwrap() < joint.evaluate(&wrong).unwrap());
    fn zero(_: &dyn Worldline) -> Result<Vec<f64>> {
        Ok(vec![0.0])
    }
    assert!(
        PhotonObservation {
            counts: vec![1],
            model: zero
        }
        .evaluate(&truth)
        .is_err()
    );
    close(
        PhotonObservation {
            counts: vec![0],
            model: zero,
        }
        .evaluate(&truth)
        .unwrap(),
        0.0,
        0.0,
    );
}

fn response(_: &ReceivedBistaticPath, _: f64) -> Result<ComplexAmplitude> {
    Ok(ComplexAmplitude {
        real: 1.0,
        imaginary: 0.0,
    })
}

#[test]
fn asteroid_echoes_arrive_together_and_can_cancel_or_reinforce() {
    let endpoint = world(Vec3::ZERO, Vec3::ZERO);
    let a = world(v(1000.0, 0.0, 0.0), Vec3::ZERO);
    let frequency = 10e6;
    let b = world(v(1000.0 + C / (4.0 * frequency), 0.0, 0.0), Vec3::ZERO);
    let c = world(v(1000.0 + C / (2.0 * frequency), 0.0, 0.0), Vec3::ZERO);
    let reception = Epoch::from_parts(4_000_000_000, 0.5).unwrap();
    let reference = received_bistatic_path(&endpoint, &a, &endpoint, reception, 1e-14)
        .unwrap()
        .emission_epoch;
    let sum = |other| {
        coherent_echo(
            CwReference {
                channel: "calibrated radar",
                frequency_hz: frequency,
                emission_epoch: reference,
            },
            &endpoint,
            &endpoint,
            reception,
            &[
                Scatterer {
                    id: "a",
                    trajectory: &a,
                    response: &response,
                },
                Scatterer {
                    id: "b",
                    trajectory: other,
                    response: &response,
                },
            ],
            1e-14,
        )
        .unwrap()
    };
    let cancelled = sum(&b);
    assert!(cancelled.amplitude.power().unwrap() < 1e-14);
    close(sum(&c).amplitude.power().unwrap(), 4.0, 1e-12);
    for contribution in cancelled.contributions {
        assert_eq!(contribution.path.reception_epoch, reception);
        assert!(contribution.path.emission_epoch < contribution.path.scattering_epoch);
        assert!(contribution.path.scattering_epoch < reception);
    }
}

#[test]
fn received_path_matches_forward_retarded_solution_for_moving_endpoints() {
    let tx = world(v(-1000.0, 0.0, 0.0), v(1.0, 2.0, 0.0));
    let target = world(v(100_000.0, 200_000.0, 0.0), v(300.0, -100.0, 0.0));
    let rx = world(v(1000.0, 500.0, 0.0), v(-20.0, 10.0, 0.0));
    let backward = received_bistatic_path(&tx, &target, &rx, epoch(10.0), 1e-14).unwrap();
    let forward = bistatic_path(&tx, &target, &rx, backward.emission_epoch, 1e-14).unwrap();
    close(
        forward
            .reception
            .epoch
            .duration_since(backward.reception_epoch),
        0.0,
        2e-14,
    );
    close(
        forward
            .scattering
            .epoch
            .duration_since(backward.scattering_epoch),
        0.0,
        2e-14,
    );
    close(backward.delay_s(), forward.total_vacuum_delay_s(), 2e-14);
}
