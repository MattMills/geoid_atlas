use geoid_atlas::calibration::*;
use geoid_atlas::observability::*;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a} != {b}");
}
fn carriers() -> Vec<Carrier> {
    [1e6, 2e6, 3e6]
        .into_iter()
        .enumerate()
        .map(|(i, f)| Carrier {
            id: format!("c{i}"),
            frequency_hz: f,
        })
        .collect()
}
fn rows() -> Vec<ClockObservation> {
    let tx = [0.1, -0.2, 0.3];
    let mut rows = vec![];
    for r in 0..2 {
        for (j, c) in carriers().iter().enumerate() {
            rows.push(ClockObservation {
                id: format!("{r}-{j}"),
                receiver: format!("r{r}"),
                carrier: c.id.clone(),
                offset_hz: tx[j]
                    + if r == 0 {
                        0.0
                    } else {
                        2.0 * c.frequency_hz * 1e-6 + 0.5
                    },
                noise_std_hz: 1.0,
            });
        }
    }
    rows
}

#[test]
fn joint_fit_recovers_relative_clocks_and_analytic_covariance() {
    let receivers = vec!["r0".into(), "r1".into()];
    let fit = fit_clock_network(&receivers, &carriers(), &rows(), "r0").unwrap();
    close(fit.receivers[1].relative_frequency_scale_ppm, 2.0, 1e-12);
    close(fit.receivers[1].relative_offset_hz, 0.5, 1e-12);
    close(fit.covariance[0][0], 1.0, 1e-12);
    close(fit.covariance[1][1], 14.0 / 3.0, 1e-12);
    close(fit.covariance[0][1], -2.0, 1e-12);
    assert_eq!(fit.identifiable_parameters, 5);
    assert_eq!(fit.degrees_of_freedom, 1);
    assert!(fit.residual_rms_hz < 1e-12);
    assert_eq!(fit.receivers[0].relative_frequency_scale_ppm, 0.0);
    for (i, expected) in [0.1, -0.2, 0.3].into_iter().enumerate() {
        close(fit.transmitters[i].offset_hz, expected, 1e-12);
    }
}

#[test]
fn changing_clock_gauge_preserves_predictions() {
    let receivers = vec!["r0".into(), "r1".into()];
    let observations = rows();
    let a = fit_clock_network(&receivers, &carriers(), &observations, "r0").unwrap();
    let b = fit_clock_network(&receivers, &carriers(), &observations, "r1").unwrap();
    close(b.receivers[0].relative_frequency_scale_ppm, -2.0, 1e-12);
    close(b.receivers[0].relative_offset_hz, -0.5, 1e-12);
    for (x, y) in a.residuals.iter().zip(&b.residuals) {
        close(x.predicted_hz, y.predicted_hz, 1e-12);
    }
}

#[test]
fn noise_weights_control_the_fit_and_covariance_units() {
    let receivers = vec!["r0".into(), "r1".into()];
    let mut observations = rows();
    observations[5].offset_hz += 3.0;
    let ordinary = fit_clock_network(&receivers, &carriers(), &observations, "r0").unwrap();
    observations[5].noise_std_hz = 100.0;
    let weighted = fit_clock_network(&receivers, &carriers(), &observations, "r0").unwrap();
    assert!(
        (weighted.receivers[1].relative_frequency_scale_ppm - 2.0).abs()
            < (ordinary.receivers[1].relative_frequency_scale_ppm - 2.0).abs() * 0.01
    );
    let base = fit_clock_network(&receivers, &carriers(), &rows(), "r0").unwrap();
    let scaled_rows: Vec<_> = rows()
        .into_iter()
        .map(|row| ClockObservation {
            noise_std_hz: 3.0,
            ..row
        })
        .collect();
    let scaled = fit_clock_network(&receivers, &carriers(), &scaled_rows, "r0").unwrap();
    close(
        scaled.receivers[1].relative_frequency_scale_ppm,
        base.receivers[1].relative_frequency_scale_ppm,
        1e-12,
    );
    close(scaled.covariance[0][0], 9.0 * base.covariance[0][0], 1e-11);
    for (i, row) in scaled.covariance.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert_eq!(*value, scaled.covariance[j][i]);
        }
    }
}

#[test]
fn rank_defects_are_errors_without_hidden_priors() {
    let receivers = vec!["r0".into(), "r1".into()];
    let mut frequencies = carriers();
    for c in &mut frequencies {
        c.frequency_hz = 1e6;
    }
    assert!(
        fit_clock_network(&receivers, &frequencies, &rows(), "r0")
            .unwrap_err()
            .to_string()
            .contains("rank deficient")
    );
    let mut observations = rows();
    for row in &mut observations {
        if row.receiver == "r1" {
            row.carrier = "c0".into();
        }
    }
    assert!(fit_clock_network(&receivers, &carriers(), &observations, "r0").is_err());
}

#[test]
fn bad_ids_and_uncertainty_are_rejected() {
    let receivers = vec!["r0".into(), "r1".into()];
    let mut observations = rows();
    observations[1].id = observations[0].id.clone();
    assert!(fit_clock_network(&receivers, &carriers(), &observations, "r0").is_err());
    let mut observations = rows();
    observations[0].noise_std_hz = 0.0;
    assert!(fit_clock_network(&receivers, &carriers(), &observations, "r0").is_err());
    let mut observations = rows();
    observations[0].carrier = "unknown".into();
    assert!(fit_clock_network(&receivers, &carriers(), &observations, "r0").is_err());
    assert!(fit_clock_network(&receivers, &carriers(), &rows(), "unknown").is_err());
}

#[test]
fn sparse_report_sized_grid_has_74_parameters_and_188_dof() {
    let frequencies = [78e3, 162e3, 540e3, 693e3, 882e3, 954e3, 1089e3, 1278e3];
    let c: Vec<_> = frequencies
        .into_iter()
        .enumerate()
        .map(|(i, f)| Carrier {
            id: format!("c{i}"),
            frequency_hz: f,
        })
        .collect();
    let r: Vec<_> = (0..34).map(|i| format!("r{i}")).collect();
    let mut observations = vec![];
    for (i, receiver) in r.iter().enumerate() {
        for (j, carrier) in c.iter().enumerate() {
            if (1..=10).contains(&i) && j == 7 {
                continue;
            }
            observations.push(ClockObservation {
                id: format!("{i}-{j}"),
                receiver: receiver.clone(),
                carrier: carrier.id.clone(),
                offset_hz: i as f64 * 0.1 * carrier.frequency_hz * 1e-6 - i as f64 * 0.02
                    + j as f64 * 0.03,
                noise_std_hz: 0.1,
            });
        }
    }
    let fit = fit_clock_network(&r, &c, &observations, "r0").unwrap();
    assert_eq!(observations.len(), 262);
    assert_eq!(fit.identifiable_parameters, 74);
    assert_eq!(fit.degrees_of_freedom, 188);
    assert!(fit.residual_rms_hz < 1e-12);
    close(fit.receivers[33].relative_frequency_scale_ppm, 3.3, 1e-10);
}

#[test]
fn kiwi_frequency_plan_and_resolution_diagnostics() {
    let frequencies = [78e3, 162e3, 540e3, 693e3, 882e3, 954e3, 1089e3, 1278e3];
    let diagnostic = clock_dispersion_separability(&frequencies, None).unwrap();
    close(diagnostic.correlation, -0.8301024807304196, 1e-12);
    close(
        diagnostic.variance_inflation.unwrap(),
        3.2161593069953343,
        1e-12,
    );
    close(
        delay_resolution(461.0, DistanceConvention::OneWayPath)
            .unwrap()
            .distance_cell_m,
        650309.0195227766,
        1e-6,
    );
    close(
        delay_resolution(1.2e6, DistanceConvention::MonostaticRange)
            .unwrap()
            .distance_cell_m,
        124.91352416666667,
        1e-9,
    );
    let confounded = clock_dispersion_separability(&[1e6, 2e6], None).unwrap();
    assert!(confounded.variance_inflation.is_none());
    assert!(clock_dispersion_separability(&[1e6, 1e6], None).is_err());
    assert!(clock_dispersion_separability(&frequencies, Some(&[1.0])).is_err());
    assert!(delay_resolution(0.0, DistanceConvention::OneWayPath).is_err());
}

#[test]
fn nominal_ambiguity_lattice_requires_phase_calibration() {
    let lattice = nominal_frequency_lattice(&[
        78_000, 162_000, 540_000, 693_000, 882_000, 954_000, 1_089_000, 1_278_000,
    ])
    .unwrap();
    assert_eq!(lattice.difference_gcd_hz, 3000);
    close(lattice.one_way_path_period_m, 99_930.81933333333, 1e-6);
    close(lattice.largest_pair_beat_path_m, 4163.784138888889, 1e-6);
    assert!(nominal_frequency_lattice(&[100, 100]).is_err());
}
