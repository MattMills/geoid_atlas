use geoid_atlas::frames::*;
use geoid_atlas::interferometry::SPEED_OF_LIGHT_M_S;
use geoid_atlas::spacetime::*;
use geoid_atlas::time::*;
use geoid_atlas::{Error, Result};
fn t(s: f64) -> Epoch {
    Epoch::new(s).unwrap()
}
fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn line(position: Vec3, velocity: Vec3) -> ConstantAcceleration {
    ConstantAcceleration {
        reference_epoch: t(0.0),
        state: State {
            position_m: position,
            velocity_m_s: velocity,
        },
        acceleration_m_s2: Vec3::ZERO,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a} != {b}");
}

#[test]
fn split_time_preserves_fine_offsets_and_requires_a_scale_model() {
    let origin = JulianDate::new(2_451_545, 40_000.0).unwrap();
    let reference = TimeReference {
        origin,
        scale: TimeScale::Tt,
    };
    let timestamp = Timestamp {
        date: origin.shifted(1e-7).unwrap(),
        scale: TimeScale::Tt,
    };
    close(
        reference.epoch_at(timestamp).unwrap().seconds(),
        1e-7,
        2e-12,
    );
    let tai = Timestamp {
        date: origin.shifted(-32.184).unwrap(),
        scale: TimeScale::Tai,
    };
    close(reference.epoch_at(tai).unwrap().seconds(), 0.0, 1e-10);
    assert!(
        Timestamp {
            date: origin,
            scale: TimeScale::Tdb
        }
        .to_scale(TimeScale::Tt)
        .is_err()
    );
    let graph = FrameGraph::with_time_reference("root", reference);
    close(graph.epoch_at(timestamp).unwrap().seconds(), 1e-7, 2e-12);
    assert!(FrameGraph::new("root").epoch_at(timestamp).is_err());
    let normalized = JulianDate::new(100, -1.0).unwrap();
    assert_eq!(normalized.day(), 99);
    assert_eq!(normalized.seconds(), 86399.0);
}

#[test]
fn fine_light_time_offsets_survive_large_epoch_origin() {
    let reference = Epoch::from_parts(900_000_000, 0.25).unwrap();
    let shifted = reference.shifted(1e-10).unwrap();
    close(shifted.duration_since(reference), 1e-10, 1e-16);
    close(
        reference.shifted(-1e-10).unwrap().duration_since(reference),
        -1e-10,
        1e-16,
    );
    let model = ConstantAcceleration {
        reference_epoch: reference,
        state: State {
            position_m: Vec3::ZERO,
            velocity_m_s: v(1e6, 0.0, 0.0),
        },
        acceleration_m_s2: Vec3::ZERO,
    };
    close(model.state_at(shifted).unwrap().position_m.x, 1e-4, 1e-10);
    let time = TimeReference {
        origin: JulianDate::new(2_451_545, 0.0).unwrap(),
        scale: TimeScale::Tdb,
    };
    let restored = time
        .epoch_at(time.timestamp_at(reference).unwrap())
        .unwrap();
    close(restored.duration_since(reference), 0.0, 1e-12);
}

#[test]
fn hermite_worldpath_recovers_constant_acceleration_and_is_contiguous() {
    let model = ConstantAcceleration {
        acceleration_m_s2: v(2.0, 0.0, -1.0),
        ..line(v(10.0, 20.0, 30.0), v(1.0, 2.0, 3.0))
    };
    let knots = [0.0, 10.0, 20.0]
        .into_iter()
        .map(|s| (t(s), model.state_at(t(s)).unwrap()))
        .collect();
    let path = HermiteWorldline::new(knots).unwrap();
    for s in [0.0, 1.0, 9.99, 10.0, 10.01, 19.0, 20.0] {
        let expected = model.state_at(t(s)).unwrap();
        let actual = path.state_at(t(s)).unwrap();
        close((actual.position_m - expected.position_m).norm(), 0.0, 1e-9);
        close(
            (actual.velocity_m_s - expected.velocity_m_s).norm(),
            0.0,
            1e-9,
        );
    }
    assert!(matches!(
        path.state_at(t(21.0)),
        Err(Error::OutsideCoverage)
    ));
    assert!(HermiteWorldline::new(vec![(t(0.0), model.state), (t(0.0), model.state)]).is_err());
}

#[test]
fn bistatic_path_uses_two_distinct_retarded_events() {
    let tx = line(Vec3::ZERO, Vec3::ZERO);
    let target = line(v(3e8, 0.0, 0.0), v(1000.0, 0.0, 0.0));
    let rx = line(Vec3::ZERO, Vec3::ZERO);
    let path = bistatic_path(&tx, &target, &rx, t(0.0), 1e-13).unwrap();
    let leg = 3e8 / (SPEED_OF_LIGHT_M_S - 1000.0);
    close(path.scattering.travel_time_s, leg, 1e-12);
    close(path.reception.travel_time_s, leg, 1e-12);
    close(path.total_vacuum_delay_s(), 2.0 * leg, 1e-12);
    assert!(
        path.reception.epoch > path.scattering.epoch && path.scattering.epoch > path.emission_epoch
    );
}

#[test]
fn rephasing_recovers_signal_and_sqrt_n_noise_gain() {
    let samples = (0..100)
        .map(|i| {
            let phase = i as f64 * 0.37;
            PhaseSample {
                real: phase.cos(),
                imaginary: phase.sin(),
                predicted_phase_rad: phase,
                noise_std: 2.0,
            }
        })
        .collect::<Vec<_>>();
    let aligned = integrate_coherent(&samples).unwrap();
    close(aligned.real, 1.0, 1e-12);
    close(aligned.imaginary, 0.0, 1e-12);
    close(aligned.noise_std, 0.2, 1e-12);
    let unaligned = samples
        .iter()
        .map(|s| PhaseSample {
            predicted_phase_rad: 0.0,
            ..*s
        })
        .collect::<Vec<_>>();
    let washed = integrate_coherent(&unaligned).unwrap();
    assert!(washed.real.hypot(washed.imaginary) < 0.1);
}

#[test]
fn coherent_trajectory_likelihood_distinguishes_wrong_worldpath() {
    let tx = [
        line(v(-1000.0, 0.0, 0.0), Vec3::ZERO),
        line(v(0.0, -1000.0, 0.0), Vec3::ZERO),
        line(v(0.0, 0.0, -1000.0), Vec3::ZERO),
    ];
    let rx = [
        line(Vec3::ZERO, Vec3::ZERO),
        line(v(100.0, 50.0, 0.0), Vec3::ZERO),
    ];
    let target = line(v(100.0, 200.0, 300.0), v(1.0, 2.0, 0.5));
    let mut observations = vec![];
    for second in [0.0, 10.0, 20.0, 30.0] {
        for (i, transmitter) in tx.iter().enumerate() {
            for (j, receiver) in rx.iter().enumerate() {
                let path = bistatic_path(transmitter, &target, receiver, t(second), 1e-13).unwrap();
                observations.push(BistaticObservation {
                    id: format!("{second}-{i}-{j}"),
                    coherent_channel: "calibrated RF".into(),
                    emission_epoch: t(second),
                    transmitter_index: i,
                    receiver_index: j,
                    frequency_hz: 30e6,
                    reference_delay_s: path.total_vacuum_delay_s(),
                    real: 1.0,
                    imaginary: 0.0,
                    noise_std: 0.1,
                });
            }
        }
    }
    let model = CoherentBistaticLikelihood::new(
        tx.iter().map(|s| s as &dyn Worldline).collect(),
        rx.iter().map(|s| s as &dyn Worldline).collect(),
        observations,
        1e-13,
    )
    .unwrap();
    close(model.evaluate(&target).unwrap(), 0.0, 1e-20);
    let wrong = line(
        target.state.position_m + v(10.0, 0.0, 0.0),
        target.state.velocity_m_s,
    );
    assert!(model.evaluate(&wrong).unwrap() > 100.0);
    let nearby = line(
        target.state.position_m + v(0.5, -0.25, 0.25),
        target.state.velocity_m_s + v(0.02, -0.01, 0.01),
    );
    let initial_cost = model.evaluate(&nearby).unwrap();
    let fit = refine_trajectory(
        nearby,
        LocalSearch {
            position_half_width_m: 1.0,
            velocity_half_width_m_s: 0.05,
            levels: 9,
            sweeps_per_level: 10,
        },
        &model,
    )
    .unwrap();
    assert!(
        fit.negative_log_likelihood < initial_cost * 0.001,
        "actual coherent-path fit did not improve enough: {} -> {}",
        initial_cost,
        fit.negative_log_likelihood
    );
}

struct PositionLikelihood {
    truth: State,
}
impl NegativeLogLikelihood for PositionLikelihood {
    fn evaluate(&self, path: &dyn Worldline) -> Result<f64> {
        let state = path.state_at(t(0.0))?;
        Ok((state.position_m - self.truth.position_m).norm().powi(2)
            + (state.velocity_m_s - self.truth.velocity_m_s)
                .norm()
                .powi(2))
    }
}
#[test]
fn bounded_local_fit_refines_six_parameters_without_double_counting() {
    let truth = State {
        position_m: v(0.5, -0.25, 0.125),
        velocity_m_s: v(0.25, 0.125, -0.5),
    };
    let initial = line(Vec3::ZERO, Vec3::ZERO);
    let fit = refine_trajectory(
        initial,
        LocalSearch {
            position_half_width_m: 1.0,
            velocity_half_width_m_s: 1.0,
            levels: 5,
            sweeps_per_level: 3,
        },
        &PositionLikelihood { truth },
    )
    .unwrap();
    close(fit.negative_log_likelihood, 0.0, 1e-12);
    assert!(fit.evaluations > 12);
    let mut joint = JointLikelihood::default();
    joint
        .add("positions", "same-survey", PositionLikelihood { truth })
        .unwrap();
    assert!(
        joint
            .add(
                "duplicate representation",
                "same-survey",
                PositionLikelihood { truth }
            )
            .is_err()
    );
    let cell = LocalCell::new(Vec3::ZERO, t(0.0), 100.0, 60.0).unwrap();
    assert!(cell.contains(v(99.0, 0.0, 0.0), t(30.0)).unwrap());
    assert!(
        !cell
            .refined()
            .unwrap()
            .contains(v(99.0, 0.0, 0.0), t(30.0))
            .unwrap()
    );
}

#[test]
fn prior_and_photon_count_likelihoods_remain_separate_models() {
    let trajectory = line(Vec3::ZERO, Vec3::ZERO);
    let prior = GaussianStatePrior {
        epoch: t(0.0),
        mean: State::stationary(v(1.0, 0.0, 0.0)),
        position_std_m: v(2.0, 2.0, 2.0),
        velocity_std_m_s: v(1.0, 1.0, 1.0),
    };
    close(prior.evaluate(&trajectory).unwrap(), 0.125, 1e-12);
    assert!(poisson_count_cost(10, 10.0).unwrap() < poisson_count_cost(10, 5.0).unwrap());
    assert_eq!(poisson_count_cost(0, 0.0).unwrap(), 0.0);
    assert!(poisson_count_cost(1, 0.0).is_err());
}
