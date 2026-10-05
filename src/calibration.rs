//! Sparse receiver/carrier frequency calibration with an explicit two-dimensional gauge.
//!
//! Model in Hz: `y_ij = epsilon_i[ppm] * f_j[Hz] * 1e-6 + b_i[Hz] + t_j[Hz]`.
//! One reference receiver has epsilon=b=0 by CONVENTION, not measurement. No
//! absolute clock, proper-time rate, phase, or GPS start-time accuracy is inferred.
//!
//! ```
//! use geoid_atlas::calibration::*;
//! let receivers: Vec<String> = vec!["reference".into(), "other".into()];
//! let carriers: Vec<_> = [1e6, 2e6, 3e6].into_iter().enumerate()
//!     .map(|(i, f)| Carrier { id: format!("c{i}"), frequency_hz: f }).collect();
//! let mut observations = Vec::new();
//! for (i, receiver) in receivers.iter().enumerate() {
//!     for carrier in &carriers {
//!         observations.push(ClockObservation {
//!             id: format!("{receiver}-{}", carrier.id), receiver: receiver.clone(),
//!             carrier: carrier.id.clone(), noise_std_hz: 0.1,
//!             offset_hz: if i == 0 { 0.0 } else { carrier.frequency_hz * 1e-6 },
//!         });
//!     }
//! }
//! let fit = fit_clock_network(&receivers, &carriers, &observations, "reference")?;
//! assert!((fit.receivers[1].relative_frequency_scale_ppm - 1.0).abs() < 1e-10);
//! # Ok::<(), geoid_atlas::Error>(())
//! ```
use crate::{Error, Result, finite};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct Carrier {
    pub id: String,
    pub frequency_hz: f64,
}
#[derive(Debug, Clone)]
pub struct ClockObservation {
    pub id: String,
    pub receiver: String,
    pub carrier: String,
    pub offset_hz: f64,
    pub noise_std_hz: f64,
}
#[derive(Debug, Clone)]
pub struct ReceiverCalibration {
    pub id: String,
    pub relative_frequency_scale_ppm: f64,
    pub relative_offset_hz: f64,
    pub std_frequency_scale_ppm: f64,
    pub std_offset_hz: f64,
}
#[derive(Debug, Clone)]
pub struct TransmitterCalibration {
    pub carrier: String,
    pub offset_hz: f64,
    pub std_offset_hz: f64,
}
#[derive(Debug, Clone)]
pub struct CalibrationResidual {
    pub observation_id: String,
    pub predicted_hz: f64,
    pub residual_hz: f64,
}
#[derive(Debug, Clone)]
pub struct ClockSolution {
    pub reference_receiver: String,
    pub receivers: Vec<ReceiverCalibration>,
    pub transmitters: Vec<TransmitterCalibration>,
    /// Covariance corresponds to parameter_names, conditional on the chosen gauge.
    /// Given noise variances are treated as known; no residual-based rescaling is applied.
    pub parameter_names: Vec<String>,
    pub covariance: Vec<Vec<f64>>,
    pub residuals: Vec<CalibrationResidual>,
    pub identifiable_parameters: usize,
    pub degrees_of_freedom: usize,
    pub weighted_chi_squared: f64,
    pub residual_rms_hz: f64,
    pub raw_rms_about_zero_hz: f64,
    pub smallest_scaled_qr_pivot: f64,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidInput(message.into())
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |total, value| total.hypot(*value))
}
struct LeastSquares {
    coefficients: Vec<f64>,
    covariance: Vec<Vec<f64>>,
    smallest_pivot: f64,
}

// Column-pivoted, twice-reorthogonalized modified Gram-Schmidt QR on scaled
// columns. Rank defects are reported; no ridge prior silently resolves a gauge.
fn weighted_least_squares(mut columns: Vec<Vec<f64>>, response: Vec<f64>) -> Result<LeastSquares> {
    let p = columns.len();
    let mut scales = vec![0.0; p];
    let mut permutation: Vec<_> = (0..p).collect();
    for (i, column) in columns.iter_mut().enumerate() {
        scales[i] = norm(column);
        finite(scales[i], "design scale")?;
        if scales[i] > 0.0 {
            for value in column {
                *value /= scales[i];
            }
        }
    }
    let mut r = vec![vec![0.0; p]; p];
    let mut qty = vec![0.0; p];
    let mut smallest = 1.0_f64;
    for k in 0..p {
        let pivot = (k..p)
            .max_by(|a, b| norm(&columns[*a]).total_cmp(&norm(&columns[*b])))
            .expect("nonempty pivot range");
        columns.swap(k, pivot);
        permutation.swap(k, pivot);
        for row in r.iter_mut().take(k) {
            row.swap(k, pivot);
        }
        let diagonal = norm(&columns[k]);
        if diagonal <= 1e-10 {
            return Err(invalid(format!(
                "clock model is rank deficient: rank {k} of {p}; check carrier diversity, sparse coverage and gauge"
            )));
        }
        smallest = smallest.min(diagonal);
        r[k][k] = diagonal;
        let q: Vec<_> = columns[k].iter().map(|value| value / diagonal).collect();
        qty[k] = dot(&q, &response);
        finite(qty[k], "QR response")?;
        for j in (k + 1)..p {
            let projection = dot(&q, &columns[j]);
            for (v, basis) in columns[j].iter_mut().zip(&q) {
                *v -= projection * basis;
            }
            let correction = dot(&q, &columns[j]);
            for (v, basis) in columns[j].iter_mut().zip(&q) {
                *v -= correction * basis;
            }
            r[k][j] = projection + correction;
        }
    }
    let mut z = vec![0.0; p];
    for i in (0..p).rev() {
        z[i] = (qty[i] - ((i + 1)..p).map(|j| r[i][j] * z[j]).sum::<f64>()) / r[i][i];
        finite(z[i], "calibration coefficient")?;
    }
    let mut coefficients = vec![0.0; p];
    for i in 0..p {
        coefficients[permutation[i]] = z[i] / scales[permutation[i]];
        finite(coefficients[permutation[i]], "calibration coefficient")?;
    }
    let inverse_columns: Vec<Vec<f64>> = (0..p)
        .map(|col| {
            let mut result = vec![0.0; p];
            for row in (0..p).rev() {
                result[row] = ((if row == col { 1.0 } else { 0.0 })
                    - ((row + 1)..p).map(|k| r[row][k] * result[k]).sum::<f64>())
                    / r[row][row];
            }
            result
        })
        .collect();
    let inverse: Vec<Vec<f64>> = (0..p)
        .map(|row| inverse_columns.iter().map(|column| column[row]).collect())
        .collect();
    let mut covariance = vec![vec![0.0; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let value =
                dot(&inverse[i], &inverse[j]) / scales[permutation[i]] / scales[permutation[j]];
            finite(value, "calibration covariance")?;
            covariance[permutation[i]][permutation[j]] = value;
            covariance[permutation[j]][permutation[i]] = value;
        }
        if covariance[permutation[i]][permutation[i]] <= 0.0 {
            return Err(invalid(
                "calibration variance underflows representable range",
            ));
        }
    }
    Ok(LeastSquares {
        coefficients,
        covariance,
        smallest_pivot: smallest,
    })
}

/// Joint frequency-scale, receiver-offset and transmitter-offset fit for one
/// stationary calibration window. Missing cells are allowed; structural rank
/// defects, unknown IDs and duplicate observation IDs are errors. Distinct rows
/// must have conditionally independent errors or be prewhitened by a joint model.
/// This is a dense local solve, not a streaming solver for an entire raw archive.
pub fn fit_clock_network(
    receivers: &[String],
    carriers: &[Carrier],
    observations: &[ClockObservation],
    reference_receiver: &str,
) -> Result<ClockSolution> {
    if receivers.is_empty() || carriers.is_empty() || observations.is_empty() {
        return Err(invalid(
            "calibration needs receivers, carriers and observations",
        ));
    }
    let mut receiver_index = HashMap::new();
    for (i, id) in receivers.iter().enumerate() {
        if id.trim().is_empty() || receiver_index.insert(id.as_str(), i).is_some() {
            return Err(invalid("receiver IDs must be nonempty and unique"));
        }
    }
    let reference = *receiver_index
        .get(reference_receiver)
        .ok_or_else(|| invalid("unknown reference receiver"))?;
    let mut carrier_index = HashMap::new();
    for (i, carrier) in carriers.iter().enumerate() {
        finite(carrier.frequency_hz, "carrier frequency")?;
        if carrier.frequency_hz <= 0.0
            || carrier.id.trim().is_empty()
            || carrier_index.insert(carrier.id.as_str(), i).is_some()
        {
            return Err(invalid("carriers need unique IDs and positive frequencies"));
        }
    }
    let mut indices = vec![None; receivers.len()];
    let mut names = vec![];
    for (i, id) in receivers.iter().enumerate() {
        if i != reference {
            indices[i] = Some(names.len());
            names.push(format!("receiver/{id}/frequency-scale-ppm"));
            names.push(format!("receiver/{id}/offset-hz"));
        }
    }
    let transmitter_start = names.len();
    for c in carriers {
        names.push(format!("transmitter/{}/offset-hz", c.id));
    }
    let parameters = names.len();
    if observations.len() < parameters {
        return Err(invalid(format!(
            "{} observations cannot identify {parameters} free parameters",
            observations.len()
        )));
    }
    parameters
        .checked_mul(observations.len())
        .ok_or_else(|| invalid("calibration design dimensions overflow"))?;
    let mut columns = vec![vec![0.0; observations.len()]; parameters];
    let mut response = vec![0.0; observations.len()];
    let mut seen = HashSet::new();
    let mut mappings = Vec::with_capacity(observations.len());
    for (row, observation) in observations.iter().enumerate() {
        if observation.id.trim().is_empty() || !seen.insert(observation.id.as_str()) {
            return Err(invalid("observation IDs must be nonempty and unique"));
        }
        let receiver = *receiver_index
            .get(observation.receiver.as_str())
            .ok_or_else(|| invalid("observation names an unknown receiver"))?;
        let carrier = *carrier_index
            .get(observation.carrier.as_str())
            .ok_or_else(|| invalid("observation names an unknown carrier"))?;
        finite(observation.offset_hz, "frequency offset")?;
        finite(observation.noise_std_hz, "frequency uncertainty")?;
        if observation.noise_std_hz <= 0.0 {
            return Err(invalid("frequency uncertainty must be positive"));
        }
        let weight = 1.0 / observation.noise_std_hz;
        finite(weight, "calibration whitening weight")?;
        if let Some(index) = indices[receiver] {
            columns[index][row] = carriers[carrier].frequency_hz * 1e-6 * weight;
            columns[index + 1][row] = weight;
        }
        columns[transmitter_start + carrier][row] = weight;
        response[row] = observation.offset_hz * weight;
        finite(response[row], "whitened frequency offset")?;
        mappings.push((receiver, carrier));
    }
    for column in &columns {
        for value in column {
            finite(*value, "calibration design entry")?;
        }
    }
    let fit = weighted_least_squares(columns, response)?;
    let receiver_solutions = receivers
        .iter()
        .enumerate()
        .map(|(i, id)| {
            if let Some(index) = indices[i] {
                ReceiverCalibration {
                    id: id.clone(),
                    relative_frequency_scale_ppm: fit.coefficients[index],
                    relative_offset_hz: fit.coefficients[index + 1],
                    std_frequency_scale_ppm: fit.covariance[index][index].sqrt(),
                    std_offset_hz: fit.covariance[index + 1][index + 1].sqrt(),
                }
            } else {
                ReceiverCalibration {
                    id: id.clone(),
                    relative_frequency_scale_ppm: 0.0,
                    relative_offset_hz: 0.0,
                    std_frequency_scale_ppm: 0.0,
                    std_offset_hz: 0.0,
                }
            }
        })
        .collect();
    let transmitters = carriers
        .iter()
        .enumerate()
        .map(|(i, c)| TransmitterCalibration {
            carrier: c.id.clone(),
            offset_hz: fit.coefficients[transmitter_start + i],
            std_offset_hz: fit.covariance[transmitter_start + i][transmitter_start + i].sqrt(),
        })
        .collect();
    let mut residuals = vec![];
    let mut chi_squared = 0.0;
    let mut residual_norm = 0.0_f64;
    let mut raw_norm = 0.0_f64;
    for (observation, (receiver, carrier)) in observations.iter().zip(mappings) {
        let predicted = indices[receiver].map_or(0.0, |index| {
            fit.coefficients[index] * carriers[carrier].frequency_hz * 1e-6
                + fit.coefficients[index + 1]
        }) + fit.coefficients[transmitter_start + carrier];
        let residual = observation.offset_hz - predicted;
        finite(predicted, "calibration prediction")?;
        finite(residual, "calibration residual")?;
        chi_squared += (residual / observation.noise_std_hz).powi(2);
        residual_norm = residual_norm.hypot(residual);
        raw_norm = raw_norm.hypot(observation.offset_hz);
        residuals.push(CalibrationResidual {
            observation_id: observation.id.clone(),
            predicted_hz: predicted,
            residual_hz: residual,
        });
    }
    finite(chi_squared, "calibration chi squared")?;
    finite(residual_norm, "calibration residual norm")?;
    finite(raw_norm, "calibration raw norm")?;
    Ok(ClockSolution {
        reference_receiver: reference_receiver.into(),
        receivers: receiver_solutions,
        transmitters,
        parameter_names: names,
        covariance: fit.covariance,
        residuals,
        identifiable_parameters: parameters,
        degrees_of_freedom: observations.len() - parameters,
        weighted_chi_squared: chi_squared,
        residual_rms_hz: residual_norm / (observations.len() as f64).sqrt(),
        raw_rms_about_zero_hz: raw_norm / (observations.len() as f64).sqrt(),
        smallest_scaled_qr_pivot: fit.smallest_pivot,
    })
}
