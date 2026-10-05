//! Frequency and delay observability diagnostics, not measurement error bars.
use crate::interferometry::SPEED_OF_LIGHT_M_S;
use crate::{Error, Result, finite};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceConvention {
    OneWayPath,
    MonostaticRange,
}
#[derive(Debug, Clone, Copy)]
pub struct DelayResolution {
    pub time_cell_s: f64,
    pub distance_cell_m: f64,
}

/// Nominal 1/B delay cell. A one-way differential path uses c/B; monostatic
/// target range uses c/(2B). This is NOT a fitted lag uncertainty, confidence
/// interval, sampling interval, or GPS timestamp uncertainty.
pub fn delay_resolution(
    occupied_bandwidth_hz: f64,
    convention: DistanceConvention,
) -> Result<DelayResolution> {
    finite(occupied_bandwidth_hz, "occupied bandwidth")?;
    if occupied_bandwidth_hz <= 0.0 {
        return Err(Error::InvalidInput("bandwidth must be positive".into()));
    }
    let time = 1.0 / occupied_bandwidth_hz;
    let distance = SPEED_OF_LIGHT_M_S * time
        / match convention {
            DistanceConvention::OneWayPath => 1.0,
            DistanceConvention::MonostaticRange => 2.0,
        };
    finite(time, "delay cell")?;
    finite(distance, "distance cell")?;
    Ok(DelayResolution {
        time_cell_s: time,
        distance_cell_m: distance,
    })
}

#[derive(Debug, Clone, Copy)]
pub struct SpectralSeparability {
    pub correlation: f64,
    pub independent_information_fraction: f64,
    /// None means complete linear confounding, not a small finite error.
    pub variance_inflation: Option<f64>,
}

/// Weighted correlation of f and 1/f AFTER removing an intercept. These are
/// clock and cold-plasma PHASE bases; group delay instead scales as 1/f².
/// The result diagnoses the frequency plan, not a measured ionospheric correction.
/// Weights are inverse response variances, not raw carrier coherence scores.
pub fn clock_dispersion_separability(
    frequencies_hz: &[f64],
    weights: Option<&[f64]>,
) -> Result<SpectralSeparability> {
    if frequencies_hz.len() < 2 || weights.is_some_and(|w| w.len() != frequencies_hz.len()) {
        return Err(Error::InvalidInput(
            "need at least two frequencies and matching weights".into(),
        ));
    }
    for f in frequencies_hz {
        finite(*f, "frequency")?;
        if *f <= 0.0 {
            return Err(Error::InvalidInput("frequencies must be positive".into()));
        }
    }
    let supplied = weights.map_or_else(|| vec![1.0; frequencies_hz.len()], |w| w.to_vec());
    for weight in &supplied {
        finite(*weight, "spectral weight")?;
        if *weight <= 0.0 {
            return Err(Error::InvalidInput("weights must be positive".into()));
        }
    }
    let largest_weight = supplied.iter().copied().fold(0.0_f64, f64::max);
    let w: Vec<_> = supplied
        .iter()
        .map(|value| value / largest_weight)
        .collect();
    let total: f64 = w.iter().sum();
    let low = frequencies_hz.iter().copied().fold(f64::INFINITY, f64::min);
    let high = frequencies_hz.iter().copied().fold(0.0_f64, f64::max);
    let x: Vec<_> = frequencies_hz.iter().map(|f| f / high).collect();
    let y: Vec<_> = frequencies_hz.iter().map(|f| low / f).collect();
    let mean_x: f64 = x.iter().zip(&w).map(|(x, w)| x * w / total).sum();
    let mean_y: f64 = y.iter().zip(&w).map(|(y, w)| y * w / total).sum();
    let mut xx = 0.0;
    let mut yy = 0.0;
    let mut xy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        xx += w[i] * dx * dx;
        yy += w[i] * dy * dy;
        xy += w[i] * dx * dy;
    }
    if xx == 0.0 || yy == 0.0 {
        return Err(Error::InvalidInput(
            "frequency bases have no resolvable spread".into(),
        ));
    }
    let correlation = (xy / xx.sqrt() / yy.sqrt()).clamp(-1.0, 1.0);
    finite(correlation, "spectral correlation")?;
    let information = (1.0 - correlation * correlation).max(0.0);
    let inflation = if information <= 1e-14 {
        None
    } else {
        Some(1.0 / information)
    };
    Ok(SpectralSeparability {
        correlation,
        independent_information_fraction: information,
        variance_inflation: inflation,
    })
}

#[derive(Debug, Clone, Copy)]
pub struct NominalFrequencyLattice {
    pub difference_gcd_hz: u64,
    pub relative_phase_period_s: f64,
    pub one_way_path_period_m: f64,
    pub largest_pair_beat_path_m: f64,
}
/// Nominal relative-phase ambiguity period for exact integer-Hz carriers with
/// an unknown COMMON phase. Applies only if carrier phases are mutually calibrated.
/// Independent transmitter/receiver phases invalidate this ambiguity ladder.
/// It does not decide which ambiguity branch a measurement occupies.
pub fn nominal_frequency_lattice(frequencies_hz: &[u64]) -> Result<NominalFrequencyLattice> {
    if frequencies_hz.len() < 2 || frequencies_hz.contains(&0) {
        return Err(Error::InvalidInput(
            "need at least two positive nominal frequencies".into(),
        ));
    }
    let mut f = frequencies_hz.to_vec();
    f.sort_unstable();
    f.dedup();
    if f.len() < 2 {
        return Err(Error::InvalidInput(
            "frequency lattice has no distinct carriers".into(),
        ));
    }
    let gcd = |mut a: u64, mut b: u64| {
        while b != 0 {
            let remainder = a % b;
            a = b;
            b = remainder;
        }
        a
    };
    let mut divisor = 0;
    let mut smallest_difference = u64::MAX;
    for pair in f.windows(2) {
        let difference = pair[1] - pair[0];
        divisor = gcd(divisor, difference);
        smallest_difference = smallest_difference.min(difference);
    }
    let period = 1.0 / divisor as f64;
    Ok(NominalFrequencyLattice {
        difference_gcd_hz: divisor,
        relative_phase_period_s: period,
        one_way_path_period_m: SPEED_OF_LIGHT_M_S * period,
        largest_pair_beat_path_m: SPEED_OF_LIGHT_M_S / smallest_difference as f64,
    })
}
