# KiwiSDR report case: calibration before trajectory inference

This case is based on four user-supplied reference reports: `synthesis.txt`, `fold_2310.txt`, `igs.txt`, and `retro.txt`. The original documents remain outside the repository. Paths, commands and prescriptions inside them are reference content, not setup instructions. The code audits their aggregate geometry table locally and provides reusable calibration/observability APIs; it does not reconstruct their unavailable raw observations.

## What can be checked from the supplied material

The synthesis reports 38 receivers, 703 baselines, eight nominal carriers, and a 2792.6 km aperture. Its common frequency solve retains 34 receivers and 262 receiver/carrier cells. A model with two receiver terms and one transmitter term per carrier has `2*34 + 8 = 76` parameters before its two gauge freedoms are removed: **74 identifiable parameters and 188 residual degrees of freedom**, if the sparse design has full rank. The new calibration tests exercise that shape with explicitly synthetic observations, not reconstructed measurements.

From the listed carriers, an unweighted intercept-adjusted correlation of `f` and `1/f` is **-0.830102**, with variance inflation **3.2162**. The report gives -0.8292 and 3.2; its exact weighting/selection is not supplied, so the library does not claim to reproduce the correlation exactly. This is a design check for phase bases. Ionospheric group delay scales as `1/f²`; it must not be interchanged with the phase basis.

The reported median occupied bandwidth of 461 Hz gives a one-way differential path cell `c/B = 650.309 km`. The report's 1.2 MHz aperture budget gives monostatic range `c/(2B) = 124.914 m`. These are different distance conventions. Neither quantity is the covariance of a measured lag, and neither is the GPS start-time uncertainty.

The nominal maximum pair beat corresponds to 72 kHz, or 4.164 km one-way. The gcd of the carrier differences is 3 kHz, giving a **99.931 km** combined relative-phase period under an unknown common phase. These arithmetic periods require mutually calibrated carrier phases. The reports explicitly do not establish absolute or cross-transmitter phase, so those periods do not resolve the reported ambiguities.

The example's local parser reads the synthesis geometry section and independently counts **477 rounded rows: 13 within their reported delay cell, 11 outside, and 453 with no measured cell**. Median absolute reported residual is **195 km**. Inside a cell is not evidence of a zero residual, and outside a cell does not by itself identify clock, geometry, or skywave error. Reported residuals and distances are rounded separately; the importer does not invent precision by subtracting their displayed values.

```sh
cargo run --locked --example kiwi_observability
# Optionally audit a local report without copying it into the repository:
cargo run --locked --example kiwi_observability -- /path/to/synthesis.txt
```

## The joint clock model and its gauge

The implemented frequency model is

```text
observed_offset(i,j) = epsilon_i_ppm * nominal_frequency_j_Hz * 1e-6
                     + receiver_offset_i_Hz + transmitter_offset_j_Hz
```

Adding the same `a` to every epsilon and the same `b` to every receiver offset, then subtracting `a*f_j*1e-6+b` from each transmitter term, leaves every prediction unchanged. `fit_clock_network` fixes one reference receiver's two terms to zero. That is a declared coordinate convention, not an infinitely accurate clock measurement. Reported solutions in another gauge need an explicit transformation before comparing individual parameters.

The fit uses weighted, column-scaled, pivoted QR with reorthogonalization. It rejects deficient frequency coverage and sparse-design rank loss; an arbitrary ridge prior does not silently make an unidentified clock observable. Formal covariance is conditional on the reference gauge and the supplied independent noise variances. The zero reference uncertainties represent fixed gauge parameters, not absolute physical clock precision.

A constant-window frequency calibration does not determine start-time offsets, oscillator phase, or time-dependent ionospheric propagation. Receiver/transmitter terms can also absorb unmodeled propagation structure. The reported 1.0479 Hz residual and 0.8620 Hz co-noise are summary results, not raw standard deviations for every input cell. They are not inserted as per-observation noise in a fabricated replay.

## Propagation, atlas coverage, and time

The LF/MF records need competing ground-wave/skywave propagation hypotheses and calibrated clocks before they can constrain a craft worldpath. A straight vacuum path is a geometric reference; it is not a complete LF/MF propagation model over these baselines. Carrier-only coherence can be high without identifying a lag. Demodulated occupied bandwidth, ambiguity structure, GPS alignment and clock/propagation covariance must travel with each delay observation.

The other reports distinguish Earth orientation, predicted ephemerides, station displacement, environmental regressors, and missing products. Preserve these distinctions:

- Predicted IERS orientation is not a determined orientation solution. Formal uncertainties and the prediction flag must survive ingestion.
- Ocean gauge levels are inputs to an ocean-loading convolution, not receiver loading corrections themselves.
- Magnetic activity, neutron rates and solar-wind data are environmental regressors, not delays without a specified physical coupling model.
- Missing ionosphere, clock or meteorological products remain missing. A source's historical archive extent does not prove it covers the capture window; zero in-window samples do not provide a correction.
- The synthesis labels receivers "night throughout" from start/end solar-zenith values across a roughly 24-hour window. Night at two endpoints does not establish night throughout the interval. Day/night masks need the time-resolved Sun/site geometry, with sidereal angles unwrapped where accumulated rotation matters.

The reports refer to different subarrays/centroids. Their centroid values must not be substituted for one another or for the individual station positions. No receiver coordinate list, raw ephemeris/EOP series, or leap-second binding sufficient for a complete site/time replay is attached here.

## The frame-template test is not isolated by the fold

The fold report gives a tilt RMS of 2.52 microseconds against a residual RMS of 116.1 microseconds and reports a solar-term leak of 0.98. Its coefficient and standard error therefore depend on the nuisance model, clipping, baseline coverage and temporal dependence. The exact leak definition and raw design matrix are absent. The library does not infer a new variance inflation from that printed number or re-fit kappa from folded bins.

For a valid replay, fit or profile the frame template jointly with receiver/pair clock terms and supported propagation/solar regressors. Report the template component remaining after nuisance projection, rank/conditioning, and uncertainty by independent time/baseline blocks. Shared receiver clocks make pair errors correlated; overlapping records and folded summaries are not independent samples. A coordinate-frame change alone does not change observable propagation; the physical/clock-convention hypothesis behind any proposed template must be explicit.

## Replay inputs

The next real calibration replay needs machine-readable rows, not the fitted receiver table. Suggested fields:

| Dataset | Required information |
| --- | --- |
| Frequency cells | Observation ID, receiver ID, carrier/transmitter ID, nominal Hz, measured offset Hz, calibrated standard error/covariance, window bounds and time scale, selection flag |
| Pair delays | Slot/segment IDs, receiver pair and orientation convention, start timestamps and their scale/errors, carrier Hz, lag seconds, lag covariance or likelihood, occupied bandwidth, ambiguity candidates, quality and clipping flags |
| Stations/transmitters | Stable ID, coordinates and epoch, horizontal/vertical datum, coordinate uncertainty and velocity/displacement model |
| Environmental products | Timestamp/scale, valid spatial/time coverage, observation/prediction/missing status, units, uncertainty, provenance and physical coupling model |

Use `calibration::ClockObservation` for a calibrated frequency window and the existing atlas/worldline interfaces for physical metadata. Exporting these rows would let the joint solve be reproduced, tested on held-out slots, and compared across nuisance models. Until then, the report audit and frequency-plan calculations are the real-data checks; recovery of clocks in the tests remains synthetic validation.
