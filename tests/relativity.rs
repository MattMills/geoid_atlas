use geoid_atlas::frames::{Epoch, State, Vec3};
use geoid_atlas::interferometry::SPEED_OF_LIGHT_M_S as C;
use geoid_atlas::relativity::{ElectromagneticField, Event, LorentzBoost};

#[test]
fn electromagnetic_boost_mixes_fields_and_preserves_both_invariants() {
    let origin = Event {
        epoch: Epoch::new(0.0).unwrap(),
        position_m: Vec3::ZERO,
    };
    let boost = LorentzBoost::new(origin, v(0.5 * C, 0.0, 0.0)).unwrap();
    let magnetic = ElectromagneticField {
        electric_v_m: Vec3::ZERO,
        magnetic_t: v(0.0, 0.0, 50e-6),
    };
    let output = boost.electromagnetic_field(magnetic).unwrap();
    let gamma = 1.0 / 0.75_f64.sqrt();
    close(output.electric_v_m.y, -gamma * 0.5 * C * 50e-6, 1e-12);
    close(output.magnetic_t.z, gamma * 50e-6, 1e-20);
    let input = ElectromagneticField {
        electric_v_m: v(1000.0, 2000.0, 3000.0),
        magnetic_t: v(10e-6, -20e-6, 30e-6),
    };
    let output = boost.electromagnetic_field(input).unwrap();
    close(
        input.electric_v_m.dot(input.magnetic_t),
        output.electric_v_m.dot(output.magnetic_t),
        1e-16,
    );
    close(
        input.magnetic_t.norm().powi(2) - input.electric_v_m.norm().powi(2) / C.powi(2),
        output.magnetic_t.norm().powi(2) - output.electric_v_m.norm().powi(2) / C.powi(2),
        1e-24,
    );
}

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn lorentz_events_preserve_interval_and_invert_with_local_epoch() {
    let origin = Event {
        epoch: Epoch::from_parts(4_000_000_000, 0.25).unwrap(),
        position_m: v(1e12, -2e12, 3e12),
    };
    let boost = LorentzBoost::new(origin, v(0.3 * C, 0.4 * C, 0.0)).unwrap();
    let event = Event {
        epoch: origin.epoch.shifted(2.0).unwrap(),
        position_m: origin.position_m + v(1e8, -2e8, 3e8),
    };
    let transformed = boost.apply(event).unwrap();
    let original_interval = C.powi(2) * 4.0 - (event.position_m - origin.position_m).norm().powi(2);
    let transformed_interval =
        C.powi(2) * transformed.epoch.seconds().powi(2) - transformed.position_m.norm().powi(2);
    close(
        original_interval / C.powi(2),
        transformed_interval / C.powi(2),
        2e-15,
    );
    let round = boost.unapply(transformed).unwrap();
    close(round.epoch.duration_since(event.epoch), 0.0, 1e-15);
    assert!((round.position_m - event.position_m).norm() < 1e-3);
}

#[test]
fn boosts_change_simultaneity_and_relativistic_velocity() {
    let origin = Event {
        epoch: Epoch::new(0.0).unwrap(),
        position_m: Vec3::ZERO,
    };
    let boost = LorentzBoost::new(origin, v(0.5 * C, 0.0, 0.0)).unwrap();
    let a = boost.apply(origin).unwrap();
    let b = boost
        .apply(Event {
            position_m: v(C, 0.0, 0.0),
            ..origin
        })
        .unwrap();
    assert!(b.epoch < a.epoch);
    let (_, state) = boost
        .apply_state(
            origin.epoch,
            State {
                position_m: Vec3::ZERO,
                velocity_m_s: v(0.75 * C, 0.0, 0.0),
            },
        )
        .unwrap();
    close(
        state.velocity_m_s.x / C,
        (0.75 - 0.5) / (1.0 - 0.75 * 0.5),
        1e-15,
    );
    let (_, rest) = boost
        .apply_state(
            origin.epoch,
            State {
                position_m: Vec3::ZERO,
                velocity_m_s: v(0.5 * C, 0.0, 0.0),
            },
        )
        .unwrap();
    assert!(rest.velocity_m_s.norm() < 1e-7);
    assert!(LorentzBoost::new(origin, v(C, 0.0, 0.0)).is_err());
}

#[test]
fn low_speed_boost_retains_small_time_shift() {
    let origin = Event {
        epoch: Epoch::from_parts(4_000_000_000, 0.0).unwrap(),
        position_m: Vec3::ZERO,
    };
    let boost = LorentzBoost::new(origin, v(100.0, 0.0, 0.0)).unwrap();
    let output = boost
        .apply(Event {
            epoch: origin.epoch,
            position_m: v(1000.0, 0.0, 0.0),
        })
        .unwrap();
    close(
        output.epoch.duration_since(Epoch::new(0.0).unwrap()),
        -100_000.0 / C / C,
        1e-16,
    );
    assert_eq!(
        LorentzBoost::new(origin, Vec3::ZERO)
            .unwrap()
            .unapply(
                LorentzBoost::new(origin, Vec3::ZERO)
                    .unwrap()
                    .apply(origin)
                    .unwrap()
            )
            .unwrap()
            .epoch,
        origin.epoch
    );
}
