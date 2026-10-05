//! Synthetic geometry only: these Kepler elements are not a real ephemeris.
use geoid_atlas::Result;
use geoid_atlas::coordinates::{Ellipsoid, Geodetic};
use geoid_atlas::frames::*;
use geoid_atlas::interferometry::far_field;
use geoid_atlas::spacetime::*;
use geoid_atlas::time::*;

fn main() -> Result<()> {
    let zero = Epoch::new(0.0)?;
    let mut graph = FrameGraph::with_time_reference(
        "synthetic inertial barycentre",
        TimeReference {
            origin: JulianDate::new(2_451_545, 0.0)?,
            scale: TimeScale::Tdb,
        },
    );
    let root = graph.root();
    let sun = graph.add("Sun-centred inertial", root, Pose::IDENTITY)?;
    let orbit = |a, e, mu, anomaly| KeplerOrbit {
        reference_epoch: zero,
        semi_major_axis_m: a,
        eccentricity: e,
        gravitational_parameter_m3_s2: mu,
        mean_anomaly_rad: anomaly,
        orbital_plane: Rotation::IDENTITY,
    };
    let earth = graph.add(
        "Earth-centred inertial",
        sun,
        orbit(149_597_870_700.0, 0.0167, 1.327_124_400_18e20, 0.0),
    )?;
    let earth_fixed = graph.add(
        "Earth-fixed",
        earth,
        LinearMotion {
            reference_epoch: zero,
            origin: State::stationary(Vec3::ZERO),
            initial_rotation: Rotation::IDENTITY,
            angular_velocity_rad_s: Vec3::new(0.0, 0.0, 7.292_115e-5)?,
        },
    )?;
    let moon = graph.add(
        "Moon-centred inertial",
        earth,
        orbit(384_400_000.0, 0.0549, 4.035_032_356e14, 1.0),
    )?;
    let satellite = graph.add(
        "Satellite array origin",
        earth,
        orbit(7_000_000.0, 0.0, 3.986_004_418e14, 0.3),
    )?;
    let target_frame = graph.add(
        "Target craft",
        earth,
        orbit(7_600_000.0, 0.01, 3.986_004_418e14, 0.9),
    )?;
    let other_star = graph.add(
        "Second star",
        root,
        Pose {
            origin: State::stationary(Vec3::new(4.0e16, 0.0, 0.0)?),
            ..Pose::IDENTITY
        },
    )?;
    let other_planet = graph.add(
        "Second star planet",
        other_star,
        orbit(1e11, 0.01, 1e20, 0.0),
    )?;
    let stationary = |p| ConstantAcceleration {
        reference_epoch: zero,
        state: State::stationary(p),
        acceleration_m_s2: Vec3::ZERO,
    };
    let transmitter_locals = (0..100)
        .map(|i| {
            let site = Geodetic::new(-5.0 + (i / 10) as f64, -5.0 + (i % 10) as f64, 10.0)?;
            let p = Ellipsoid::WGS84.to_cartesian(site);
            Ok(stationary(Vec3::new(p.x, p.y, p.z)?))
        })
        .collect::<Result<Vec<_>>>()?;
    let target_local = stationary(Vec3::ZERO);
    let target = FramedWorldline {
        graph: &graph,
        frame: target_frame,
        local: &target_local,
    };
    let receivers_local = [
        stationary(Vec3::ZERO),
        stationary(Vec3::new(100.0, 0.0, 0.0)?),
    ];
    let receivers = receivers_local
        .iter()
        .map(|local| FramedWorldline {
            graph: &graph,
            frame: satellite,
            local,
        })
        .collect::<Vec<_>>();
    let mut paths = 0;
    let mut min_delay = f64::INFINITY;
    let mut max_delay = 0.0_f64;
    for minute in 0..=60 {
        let epoch = Epoch::new(f64::from(minute) * 60.0)?;
        for local in &transmitter_locals {
            let transmitter = FramedWorldline {
                graph: &graph,
                frame: earth_fixed,
                local,
            };
            for receiver in &receivers {
                let path = bistatic_path(&transmitter, &target, receiver, epoch, 1e-10)?;
                let delay = path.total_vacuum_delay_s();
                min_delay = min_delay.min(delay);
                max_delay = max_delay.max(delay);
                paths += 1;
            }
        }
    }
    assert_eq!(paths, 12_200);
    assert!(max_delay > min_delay && min_delay > 0.0);
    println!("100 Earth transmitters × 2 orbital receivers × 61 epochs: {paths} geometric paths");
    println!("Bistatic vacuum delay range: {min_delay:.9} .. {max_delay:.9} s");
    for second in [0.0, 3600.0] {
        let epoch = Epoch::new(second)?;
        for frame in [earth, earth_fixed, moon, satellite, other_planet] {
            let state = graph.transform(State::stationary(Vec3::ZERO), sun, frame, epoch)?;
            println!(
                "t={second:.0}s Sun in {}: ({:.6e}, {:.6e}, {:.6e}) m",
                graph.name(frame)?,
                state.position_m.x,
                state.position_m.y,
                state.position_m.z
            );
        }
        let a = receivers[0].state_at(epoch)?;
        let b = receivers[1].state_at(epoch)?;
        let sun_position = graph.pose_in_root(sun, epoch)?.origin.position_m;
        let delay = far_field(a, b, sun_position - a.position_m)?;
        println!(
            "Sun far-field array delay at t={second:.0}s: {:.12e} s",
            delay.seconds
        );
    }
    println!(
        "Synthetic geometry only; no measured signals, occultation, or precision ephemerides included."
    );
    Ok(())
}
