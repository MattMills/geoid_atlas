//! Synthetic hour-long timing track. Six independent timing channels constrain
//! one local inertial target; later observations refine the earlier history.
use geoid_atlas::Result;
use geoid_atlas::estimation::*;
use geoid_atlas::frames::*;
use geoid_atlas::spacetime::*;

// Reproducible synthetic Gaussian samples. This is not an acquisition model.
struct Noise(u64);
impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 11) as f64 + 0.5) / 9_007_199_254_740_992.0
    }
    fn gaussian(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (2.0 * std::f64::consts::PI * self.uniform()).cos()
    }
}
fn main() -> Result<()> {
    let start = Epoch::from_parts(900_000_000, 0.25)?;
    let stationary = |p| ConstantAcceleration {
        reference_epoch: start,
        state: State::stationary(p),
        acceleration_m_s2: Vec3::ZERO,
    };
    let transmitters = [
        stationary(Vec3::new(-20_000.0, 0.0, 0.0)?),
        stationary(Vec3::new(0.0, -20_000.0, 0.0)?),
        stationary(Vec3::new(0.0, 0.0, -20_000.0)?),
    ];
    let receivers = [
        stationary(Vec3::ZERO),
        stationary(Vec3::new(10_000.0, 5_000.0, 0.0)?),
    ];
    let truth = ConstantAcceleration {
        reference_epoch: start,
        state: State {
            position_m: Vec3::new(100.0, 200.0, 300.0)?,
            velocity_m_s: Vec3::new(1.0, 2.0, 0.5)?,
        },
        acceleration_m_s2: Vec3::ZERO,
    };
    let initial = TrackState::new(
        start,
        State {
            position_m: truth.state.position_m + Vec3::new(10.0, -7.0, 5.0)?,
            velocity_m_s: truth.state.velocity_m_s + Vec3::new(0.1, -0.2, 0.1)?,
        },
        Covariance6::isotropic(20.0, 1.0)?,
    )?;
    let mut tracker = RecursiveTracker::new(initial, Vec3::new(1e-8, 1e-8, 1e-8)?)?;
    let mut noise = Noise(0x73b1_409e_ef24_a178);
    let mut history = vec![];
    let mut predictions = vec![];
    let mut rejected = 0;
    for minute in 0..=60 {
        let epoch = start.shifted(f64::from(minute) * 60.0)?;
        for (i, tx) in transmitters.iter().enumerate() {
            for (j, rx) in receivers.iter().enumerate() {
                let delay = bistatic_path(tx, &truth, rx, epoch, 1e-16)?.total_vacuum_delay_s();
                let outlier = if minute == 30 && i == 0 && j == 0 {
                    1e-4
                } else {
                    0.0
                };
                let observation = BistaticDelayMeasurement {
                    transmitter: tx,
                    receiver: rx,
                    emission_epoch: epoch,
                    observed_delay_s: delay + 1e-8 * noise.gaussian() + outlier,
                    noise_std_s: 1e-8,
                    known_delay_offset_s: 0.0,
                    position_step_m: 1.0,
                    velocity_step_m_s: 1.0,
                    tolerance_s: 1e-16,
                };
                let record = tracker.assimilate(
                    format!("{minute}-{i}-{j}"),
                    epoch,
                    &observation,
                    Some(25.0),
                )?;
                if minute > 0 && i == 0 && j == 0 {
                    predictions.push(record.predicted);
                }
                if !record.diagnostics.accepted {
                    rejected += 1;
                }
            }
        }
        history.push(tracker.current());
    }
    let smoothed = smooth_history(&history, &predictions)?;
    for (filtered, smooth) in history.iter().zip(&smoothed) {
        for (before, after) in filtered
            .covariance
            .standard_deviations()
            .iter()
            .zip(smooth.covariance.standard_deviations())
        {
            assert!(after <= before + 1e-8);
        }
    }
    let rms = |states: &[TrackState]| -> Result<f64> {
        let mut squared = 0.0;
        for state in states {
            squared += (state.mean.position_m - truth.state_at(state.epoch)?.position_m)
                .norm()
                .powi(2);
        }
        Ok((squared / states.len() as f64).sqrt())
    };
    let filtered_rms = rms(&history)?;
    let smoothed_rms = rms(&smoothed)?;
    assert!(smoothed_rms < filtered_rms);
    assert!(smoothed_rms < 2.0);
    assert_eq!(rejected, 1);
    let path = mean_worldpath(&smoothed)?;
    let midpoint = path.state_at(start.shifted(1800.5)?)?;
    println!("366 timing observations over one hour; {rejected} injected outlier rejected");
    println!(
        "Synthetic 3D position RMS: filtered {filtered_rms:.3} m; smoothed {smoothed_rms:.3} m"
    );
    println!(
        "Earlier position uncertainty at first epoch: filtered {:?}; smoothed {:?}",
        &history[0].covariance.standard_deviations()[..3],
        &smoothed[0].covariance.standard_deviations()[..3]
    );
    println!(
        "Contiguous interpolated mean at t=1800.5 s: {:?}",
        midpoint.position_m
    );
    println!(
        "Synthetic local inertial scene; no real acquisition, calibration, or ephemeris data."
    );
    Ok(())
}
