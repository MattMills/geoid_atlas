//! Synthetic polar calibration payload: GNSS/thermometer/pressure, aircraft
//! bearing, radar reflection and tagged photon counts share a moving worldpath.
//! Environmental profiles are conditioned inputs, not simultaneously estimated.
use geoid_atlas::Result;
use geoid_atlas::coordinates::{Ecef, Ellipsoid, Geodetic};
use geoid_atlas::echo::*;
use geoid_atlas::environment::*;
use geoid_atlas::frames::*;
use geoid_atlas::fusion::*;
use geoid_atlas::reference::tangent_pose;
use geoid_atlas::rf::ionospheric_delay;
use geoid_atlas::spacetime::*;

struct SyntheticPlasma;
impl ElectronDensity for SyntheticPlasma {
    fn per_m3_at(&self, p: Vec3, _: Epoch) -> Result<f64> {
        let approximate_height = p.norm() - 6_371_000.0;
        Ok(5e11 * (-0.5 * ((approximate_height - 350_000.0) / 100_000.0).powi(2)).exp())
    }
}
fn local(position: Vec3) -> ConstantAcceleration {
    ConstantAcceleration {
        reference_epoch: Epoch::new(0.0).unwrap(),
        state: State::stationary(position),
        acceleration_m_s2: Vec3::ZERO,
    }
}
fn response(_: &ReceivedBistaticPath, _: f64) -> Result<ComplexAmplitude> {
    // Calibrated unit response, deliberately no invented target radar cross-section.
    Ok(ComplexAmplitude {
        real: 1.0,
        imaginary: 0.0,
    })
}

fn main() -> Result<()> {
    let mut graph = FrameGraph::new("common inertial demo");
    let body = graph.add(
        "rotating Earth demo",
        graph.root(),
        LinearMotion {
            reference_epoch: Epoch::new(0.0)?,
            origin: State::stationary(Vec3::ZERO),
            initial_rotation: Rotation::IDENTITY,
            angular_velocity_rad_s: Vec3::new(0.0, 0.0, 7.292_115e-5)?,
        },
    )?;
    let site = Geodetic::new(80.0, -40.0, 0.0)?;
    let enu = graph.add(
        "Arctic station ENU",
        body,
        tangent_pose(Ellipsoid::WGS84, site)?,
    )?;
    let station_local = local(Vec3::ZERO);
    let aircraft_local = local(Vec3::new(-5000.0, 0.0, 3000.0)?);
    let truth_local = local(Vec3::new(1000.0, 500.0, 10_000.0)?);
    let wrong_local = local(Vec3::new(1100.0, 500.0, 10_300.0)?);
    let station = FramedWorldline {
        graph: &graph,
        frame: enu,
        local: &station_local,
    };
    let aircraft = FramedWorldline {
        graph: &graph,
        frame: enu,
        local: &aircraft_local,
    };
    let truth = FramedWorldline {
        graph: &graph,
        frame: enu,
        local: &truth_local,
    };
    let wrong = FramedWorldline {
        graph: &graph,
        frame: enu,
        local: &wrong_local,
    };
    let time = Epoch::new(60.0)?;
    let dipole = DipoleField {
        centre_m: Vec3::ZERO,
        moment_a_m2: Vec3::new(0.0, 0.0, 8e22)?,
        minimum_radius_m: 6e6,
    };
    let field = FramedMagneticField {
        graph: &graph,
        frame: body,
        field: &dipole,
    };
    let station_position = station.state_at(time)?.position_m;
    let axes = field_aligned_basis(
        field.teslas_at(station_position, time)?,
        Vec3::new(1.0, 0.0, 0.0)?,
    )?;
    let b_local = axes
        .inverse()
        .apply(field.teslas_at(station_position, time)?);
    assert!(b_local.x.abs() < 1e-15 && b_local.y.abs() < 1e-15);
    let outward = station_position.normalized()?;
    let gps_position = station_position + outward * 20_200_000.0;
    let gps_light_time = 20_200_000.0 / geoid_atlas::interferometry::SPEED_OF_LIGHT_M_S;
    let plasma = integrate_plasma(
        &SyntheticPlasma,
        &field,
        station_position,
        gps_position,
        time,
        time.shifted(gps_light_time)?,
        1024,
    )?;
    let gps_delay = ionospheric_delay(plasma.electron_column_per_m2, 1_575_420_000.0)?;

    // A calibrated synthetic payload records weather + GNSS height. Balloon and
    // ionosonde measurements would condition this profile with their own coverage.
    let weather_gnss = |hypothesis: &dyn Worldline| -> Result<Vec<f64>> {
        let p = graph
            .transform(hypothesis.state_at(time)?, graph.root(), body, time)?
            .position_m;
        let h = Ellipsoid::WGS84
            .to_geodetic(Ecef {
                x: p.x,
                y: p.y,
                z: p.z,
            })?
            .height_m();
        Ok(vec![
            288.15 - 0.0065 * h,
            101_325.0 * (-h / 8500.0).exp(),
            h,
        ])
    };
    let observed = weather_gnss(&truth)?;
    let gaussian = GaussianObservation::new(
        observed,
        vec![
            vec![0.25, 10.0, 0.0],
            vec![10.0, 10_000.0, 0.0],
            vec![0.0, 0.0, 25.0],
        ],
        weather_gnss,
    )?;
    let aircraft_position = aircraft.state_at(time)?.position_m;
    let emitted = emission_epoch(&truth, aircraft_position, time, 1e-13)?;
    let bearing = BearingObservation {
        sensor: &aircraft,
        reception_epoch: time,
        observed_direction: (truth.state_at(emitted)?.position_m - aircraft_position)
            .normalized()?,
        angular_std_rad: 1e-4,
        light_time_tolerance_s: 1e-13,
    };
    let path = received_bistatic_path(&station, &truth, &station, time, 1e-13)?;
    let rf_model = |hypothesis: &dyn Worldline| -> Result<Vec<f64>> {
        let echo = coherent_echo(
            CwReference {
                channel: "tagged payload RF",
                frequency_hz: 200e6,
                emission_epoch: path.emission_epoch,
            },
            &station,
            &station,
            time,
            &[Scatterer {
                id: "payload",
                trajectory: hypothesis,
                response: &response,
            }],
            1e-13,
        )?;
        Ok(vec![echo.amplitude.real, echo.amplitude.imaginary])
    };
    let rf = GaussianObservation::new(
        rf_model(&truth)?,
        vec![vec![0.01, 0.0], vec![0.0, 0.01]],
        rf_model,
    )?;
    // Separate photon likelihood; calibrated test source/response, not X-ray/RF phase locking.
    let reference_distance = (truth
        .state_at(emission_epoch(&truth, station_position, time, 1e-13)?)?
        .position_m
        - station_position)
        .norm();
    let photons = PhotonObservation {
        counts: vec![20],
        model: |hypothesis: &dyn Worldline| -> Result<Vec<f64>> {
            let emitted = emission_epoch(hypothesis, station_position, time, 1e-13)?;
            let distance = (hypothesis.state_at(emitted)?.position_m - station_position).norm();
            Ok(vec![20.0 * (reference_distance / distance).powi(2)])
        },
    };
    let mut joint = JointLikelihood::default();
    joint.add(
        "payload weather/GNSS",
        "joint payload calibration",
        gaussian,
    )?;
    joint.add("aircraft image", "camera", bearing)?;
    joint.add("RF reflection", "radar calibration", rf)?;
    joint.add("X-ray counts", "photon counter", photons)?;
    let correct = joint.evaluate(&truth)?;
    let incorrect = joint.evaluate(&wrong)?;
    assert!(incorrect > correct + 100.0);
    println!(
        "Synthetic Arctic field magnitude: {:.3} microtesla",
        b_local.z * 1e6
    );
    println!(
        "Synthetic GNSS column: {:.3} TECU; group delay {:.3} ns; Faraday {:.6} rad",
        plasma.electron_column_per_m2 / 1e16,
        gps_delay.group_delay_s * 1e9,
        plasma.faraday_rotation_rad(1_575_420_000.0)?
    );
    println!("Joint fixed-environment cost: correct {correct:.3}, displaced {incorrect:.3}");
    Ok(())
}
