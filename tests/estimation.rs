use geoid_atlas::Result;
use geoid_atlas::estimation::*;
use geoid_atlas::frames::*;
use geoid_atlas::spacetime::*;

fn epoch(seconds: f64) -> Epoch {
    Epoch::new(seconds).unwrap()
}
fn vector(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} differs from {b}");
}
fn track() -> TrackState {
    TrackState::new(
        epoch(0.0),
        State::stationary(Vec3::ZERO),
        Covariance6::isotropic(1.0, 1.0).unwrap(),
    )
    .unwrap()
}
fn x_observation(observed: f64, noise_std: f64) -> LinearMeasurement {
    LinearMeasurement {
        coefficients: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        observed,
        noise_std,
    }
}

#[test]
fn covariance_checks_full_psd_and_mixed_unit_scales() {
    assert!(Covariance6::diagonal([-1.0, 1.0, 1.0, 1.0, 1.0, 1.0]).is_err());
    let mut invalid = Covariance6::isotropic(1.0, 1.0).unwrap().matrix();
    invalid[0][1] = 2.0;
    invalid[1][0] = 2.0;
    assert!(Covariance6::new(invalid).is_err());
    invalid[0][1] = 0.2;
    invalid[1][0] = 0.1;
    assert!(Covariance6::new(invalid).is_err());
    let covariance = Covariance6::diagonal([1e20, 1e20, 1e20, 1e-20, 1e-20, 1e-20]).unwrap();
    close(covariance.standard_deviations()[3], 1e-10, 1e-20);
    let mut singular = [[0.0; 6]; 6];
    singular[0][0] = 1.0;
    singular[1][1] = 1.0;
    singular[0][1] = 1.0;
    singular[1][0] = 1.0;
    assert!(Covariance6::new(singular).is_ok());
    assert!(Covariance6::isotropic(0.0, 0.0).is_ok());
}

#[test]
fn continuous_white_acceleration_prediction_has_correct_cross_terms() {
    let state = TrackState::new(
        epoch(0.0),
        State {
            position_m: vector(1.0, 2.0, 3.0),
            velocity_m_s: vector(4.0, 5.0, 6.0),
        },
        Covariance6::isotropic(2.0, 3.0).unwrap(),
    )
    .unwrap();
    let predicted = state.predict(epoch(2.0), vector(3.0, 0.0, 0.0)).unwrap();
    let p = predicted.covariance.matrix();
    close(predicted.mean.position_m.x, 9.0, 1e-12);
    close(predicted.mean.velocity_m_s.x, 4.0, 1e-12);
    close(p[0][0], 4.0 + 4.0 * 9.0 + 3.0 * 8.0 / 3.0, 1e-12);
    close(p[0][3], 2.0 * 9.0 + 3.0 * 4.0 / 2.0, 1e-12);
    close(p[3][3], 9.0 + 3.0 * 2.0, 1e-12);
    close(p[1][1], 40.0, 1e-12);
    close(p[1][4], 18.0, 1e-12);
    assert!(state.predict(epoch(-1.0), Vec3::ZERO).is_err());
    assert!(state.predict(epoch(1.0), vector(-1.0, 0.0, 0.0)).is_err());
}

#[test]
fn kalman_update_matches_closed_form_and_gates_outliers() {
    let predicted = track().predict(epoch(1.0), Vec3::ZERO).unwrap();
    let measurement = x_observation(1.0, 1.0).linearize(&predicted).unwrap();
    let (filtered, diagnostics) = predicted.update(measurement, None).unwrap();
    let p = filtered.covariance.matrix();
    close(filtered.mean.position_m.x, 2.0 / 3.0, 1e-12);
    close(filtered.mean.velocity_m_s.x, 1.0 / 3.0, 1e-12);
    close(p[0][0], 2.0 / 3.0, 1e-12);
    close(p[0][3], 1.0 / 3.0, 1e-12);
    close(p[3][3], 2.0 / 3.0, 1e-12);
    close(diagnostics.innovation_variance, 3.0, 1e-12);
    close(diagnostics.normalized_innovation_squared, 1.0 / 3.0, 1e-12);
    let bad = x_observation(1000.0, 1.0).linearize(&predicted).unwrap();
    let (rejected, diagnostics) = predicted.update(bad, Some(9.0)).unwrap();
    assert!(!diagnostics.accepted);
    assert_eq!(rejected, predicted);
    assert!(predicted.update(measurement, Some(-1.0)).is_err());
}

#[test]
fn frame_covariance_includes_spin_and_round_trips_with_shared_distant_parent() {
    let mut graph = FrameGraph::new("inertial");
    let root = graph.root();
    let distant = graph
        .add(
            "distant star",
            root,
            Pose {
                origin: State::stationary(vector(1e20, 0.0, 0.0)),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let rotating = graph
        .add(
            "rotating",
            distant,
            Pose {
                angular_velocity_rad_s: vector(0.0, 0.0, 2.0),
                ..Pose::IDENTITY
            },
        )
        .unwrap();
    let original = TrackState::new(
        epoch(0.0),
        State::stationary(vector(0.001, 100.0, 0.0)),
        Covariance6::isotropic(2.0, 3.0).unwrap(),
    )
    .unwrap();
    let transformed = original.transform(&graph, rotating, distant).unwrap();
    let p = transformed.covariance.matrix();
    close(p[3][3], 25.0, 1e-12);
    close(p[4][4], 25.0, 1e-12);
    close(p[5][5], 9.0, 1e-12);
    close(p[0][4], 8.0, 1e-12);
    close(p[1][3], -8.0, 1e-12);
    let back = transformed.transform(&graph, distant, rotating).unwrap();
    close(
        (back.mean.position_m - original.mean.position_m).norm(),
        0.0,
        1e-12,
    );
    close(
        (back.mean.velocity_m_s - original.mean.velocity_m_s).norm(),
        0.0,
        1e-12,
    );
    for (i, row) in back.covariance.matrix().iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            close(*value, original.covariance.matrix()[i][j], 1e-12);
        }
    }
}

#[test]
fn pose_inverse_matches_unapply_for_nested_motion() {
    let pose = Pose {
        origin: State {
            position_m: vector(100.0, -20.0, 10.0),
            velocity_m_s: vector(2.0, 3.0, 4.0),
        },
        rotation: Rotation::axis_angle(vector(1.0, 2.0, 3.0), 0.7).unwrap(),
        angular_velocity_rad_s: vector(0.01, 0.02, 0.03),
    };
    let state = State {
        position_m: vector(10.0, 30.0, 40.0),
        velocity_m_s: vector(1.0, 2.0, 3.0),
    };
    let direct = pose.unapply(state);
    let inverse = pose.inverse().apply(state);
    close((direct.position_m - inverse.position_m).norm(), 0.0, 1e-12);
    close(
        (direct.velocity_m_s - inverse.velocity_m_s).norm(),
        0.0,
        1e-12,
    );
}

#[test]
fn sibling_frame_pose_and_state_jacobian_match_direct_transforms() {
    let mut graph = FrameGraph::new("root");
    let root = graph.root();
    let pose = |angle: f64, position: Vec3, spin: Vec3| Pose {
        origin: State {
            position_m: position,
            velocity_m_s: vector(1.0, 2.0, 3.0),
        },
        rotation: Rotation::axis_angle(vector(1.0, 2.0, 3.0), angle).unwrap(),
        angular_velocity_rad_s: spin,
    };
    let a = graph
        .add(
            "a",
            root,
            pose(0.3, vector(100.0, 20.0, 30.0), vector(0.01, 0.02, 0.03)),
        )
        .unwrap();
    let source = graph
        .add(
            "source",
            a,
            pose(0.7, vector(5.0, 6.0, 7.0), vector(0.03, 0.01, 0.02)),
        )
        .unwrap();
    let b = graph
        .add(
            "b",
            root,
            pose(-0.4, vector(-20.0, 40.0, 10.0), vector(-0.01, 0.04, 0.03)),
        )
        .unwrap();
    let target = graph
        .add(
            "target",
            b,
            pose(-0.9, vector(2.0, 3.0, 4.0), vector(0.01, -0.03, 0.02)),
        )
        .unwrap();
    let state = State {
        position_m: vector(20.0, 10.0, 30.0),
        velocity_m_s: vector(4.0, 5.0, 6.0),
    };
    let direct = graph
        .pose_in_root(target, epoch(0.0))
        .unwrap()
        .unapply(graph.pose_in_root(source, epoch(0.0)).unwrap().apply(state));
    let relative = graph.relative_pose(source, target, epoch(0.0)).unwrap();
    let actual = relative.apply(state);
    close((direct.position_m - actual.position_m).norm(), 0.0, 1e-12);
    close(
        (direct.velocity_m_s - actual.velocity_m_s).norm(),
        0.0,
        1e-12,
    );
    let jacobian = pose_jacobian(relative).unwrap();
    for (col, _) in jacobian[0].iter().enumerate() {
        let mut plus = state;
        let mut minus = state;
        let step = 1e-4;
        let offset = match col % 3 {
            0 => vector(step, 0.0, 0.0),
            1 => vector(0.0, step, 0.0),
            _ => vector(0.0, 0.0, step),
        };
        if col < 3 {
            plus.position_m = plus.position_m + offset;
            minus.position_m = minus.position_m - offset;
        } else {
            plus.velocity_m_s = plus.velocity_m_s + offset;
            minus.velocity_m_s = minus.velocity_m_s - offset;
        }
        let high = relative.apply(plus);
        let low = relative.apply(minus);
        let dp = (high.position_m - low.position_m) * (1.0 / (2.0 * step));
        let dv = (high.velocity_m_s - low.velocity_m_s) * (1.0 / (2.0 * step));
        for (row, derivative) in [dp.x, dp.y, dp.z, dv.x, dv.y, dv.z].into_iter().enumerate() {
            close(jacobian[row][col], derivative, 2e-9);
        }
    }
}

#[test]
fn tracker_rejects_duplicates_old_epochs_and_failed_updates_transactionally() {
    let initial = track();
    let mut tracker = RecursiveTracker::new(initial, Vec3::ZERO).unwrap();
    let measurement = x_observation(1.0, 1.0);
    tracker
        .assimilate("one", epoch(1.0), &measurement, None)
        .unwrap();
    let current = tracker.current();
    assert!(
        tracker
            .assimilate("one", epoch(2.0), &measurement, None)
            .is_err()
    );
    assert_eq!(tracker.current(), current);
    assert!(
        tracker
            .assimilate("old", epoch(0.0), &measurement, None)
            .is_err()
    );
    assert_eq!(tracker.current(), current);
    assert!(
        tracker
            .assimilate("retry", epoch(2.0), &x_observation(f64::NAN, 1.0), None)
            .is_err()
    );
    assert_eq!(tracker.current(), current);
    assert!(
        tracker
            .assimilate("retry", epoch(2.0), &measurement, None)
            .is_ok()
    );
    let rejected = tracker
        .assimilate(
            "outlier",
            epoch(3.0),
            &x_observation(1000.0, 1.0),
            Some(9.0),
        )
        .unwrap();
    assert!(!rejected.diagnostics.accepted);
    assert!(
        tracker
            .assimilate("outlier", epoch(3.0), &measurement, None)
            .is_err()
    );
}

#[test]
fn historical_smoothing_matches_two_epoch_analytic_posterior() {
    let initial = track();
    let predicted = initial.predict(epoch(1.0), Vec3::ZERO).unwrap();
    let filtered = predicted
        .update(x_observation(1.0, 1.0).linearize(&predicted).unwrap(), None)
        .unwrap()
        .0;
    let history = smooth_history(&[initial, filtered], &[predicted]).unwrap();
    let p = history[0].covariance.matrix();
    close(history[0].mean.position_m.x, 1.0 / 3.0, 1e-12);
    close(history[0].mean.velocity_m_s.x, 1.0 / 3.0, 1e-12);
    close(p[0][0], 2.0 / 3.0, 1e-12);
    close(p[0][3], -1.0 / 3.0, 1e-12);
    close(p[3][3], 2.0 / 3.0, 1e-12);
    assert_eq!(history[1], filtered);
    let path = mean_worldpath(&history).unwrap();
    let middle = path.state_at(epoch(0.5)).unwrap();
    close(middle.position_m.x, 0.5, 1e-12);
    close(middle.velocity_m_s.x, 1.0 / 3.0, 1e-12);
}

#[test]
fn smoother_checks_history_order_transition_and_singular_predictions() {
    let initial = track();
    let predicted = initial.predict(epoch(1.0), Vec3::ZERO).unwrap();
    assert!(smooth_history(&[], &[]).is_err());
    assert!(smooth_history(&[initial, predicted], &[]).is_err());
    assert!(smooth_history(&[initial, initial], &[initial]).is_err());
    let mut wrong = predicted;
    wrong.mean.position_m.x = 1.0;
    assert!(smooth_history(&[initial, predicted], &[wrong]).is_err());
    let singular = TrackState::new(
        epoch(0.0),
        State::stationary(Vec3::ZERO),
        Covariance6::isotropic(0.0, 0.0).unwrap(),
    )
    .unwrap();
    let next = singular.predict(epoch(1.0), Vec3::ZERO).unwrap();
    assert!(smooth_history(&[singular, next], &[next]).is_err());
}

#[test]
fn retarded_timing_jacobian_matches_stationary_target_geometry() {
    let tx = ConstantAcceleration {
        reference_epoch: epoch(0.0),
        state: State::stationary(Vec3::ZERO),
        acceleration_m_s2: Vec3::ZERO,
    };
    let target = TrackState::new(
        epoch(0.0),
        State::stationary(vector(1e6, 0.0, 0.0)),
        Covariance6::isotropic(100.0, 10.0).unwrap(),
    )
    .unwrap();
    let measurement = BistaticDelayMeasurement {
        transmitter: &tx,
        receiver: &tx,
        emission_epoch: epoch(0.0),
        observed_delay_s: 2e6 / 299_792_458.0,
        noise_std_s: 1e-8,
        known_delay_offset_s: 0.0,
        position_step_m: 1.0,
        velocity_step_m_s: 1.0,
        tolerance_s: 1e-16,
    };
    let linear = measurement.linearize(&target).unwrap();
    close(linear.predicted, linear.observed, 1e-15);
    close(linear.jacobian[0], 2.0 / 299_792_458.0, 1e-17);
    close(linear.jacobian[3], 2e6 / 299_792_458.0_f64.powi(2), 1e-17);
    close(linear.jacobian[1], 0.0, 1e-17);
    close(linear.jacobian[4], 0.0, 1e-17);
    let mut shifted = target;
    shifted.mean.position_m.x += 10.0;
    let observed = measurement.linearize(&shifted).unwrap();
    let filtered = shifted.update(observed, None).unwrap().0;
    assert!((filtered.mean.position_m.x - 1e6).abs() < 0.03);
    assert!(
        measurement
            .linearize(&TrackState {
                epoch: epoch(1.0),
                ..target
            })
            .is_err()
    );
}

#[test]
fn position_update_propagates_existing_cross_axis_correlation() {
    let mut covariance = Covariance6::isotropic(1.0, 1.0).unwrap().matrix();
    covariance[0][1] = 0.5;
    covariance[1][0] = 0.5;
    let prior = TrackState::new(
        epoch(0.0),
        State::stationary(Vec3::ZERO),
        Covariance6::new(covariance).unwrap(),
    )
    .unwrap();
    let measurement = x_observation(2.0, 1.0).linearize(&prior).unwrap();
    let posterior = prior.update(measurement, None).unwrap().0;
    close(posterior.mean.position_m.x, 1.0, 1e-12);
    close(posterior.mean.position_m.y, 0.5, 1e-12);
}

#[test]
fn zero_information_measurement_preserves_state() -> Result<()> {
    let initial = track();
    let model = LinearMeasurement {
        coefficients: [0.0; 6],
        observed: 2.0,
        noise_std: 1.0,
    };
    let (posterior, diagnostics) = initial.update(model.linearize(&initial)?, None)?;
    assert_eq!(initial, posterior);
    close(diagnostics.normalized_innovation_squared, 4.0, 1e-12);
    Ok(())
}
