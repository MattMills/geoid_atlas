# Extreme coordinate and sensor-fusion scenarios

The central object is a physical hypothesis in space and coordinate time. Coordinates are representations of its events; measurements constrain it through calibrated forward models. Re-expressing one measurement in six frames does not create six observations. These scenarios describe executable building blocks and the inputs needed for a credible solve, rather than claiming bundled global data or universal automatic conversion.

## Scenarios and implementation paths

| Scenario pushed to its extreme | Required inputs | Library implementation and boundary |
| --- | --- | --- |
| A buried antenna surveyed in US survey feet, a WGS84 GNSS receiver, and an orthometric DEM on another terrestrial datum | Horizontal reference realization and survey epoch; original units; both vertical models; actual datum shift/grid; antenna velocity and covariance | `units`, `GeodeticReference`, `Helmert`, raster/atlas. Convert H to h using the source geoid, transform Cartesian coordinates, invert the target ellipsoid, then use the target geoid. A Helmert cannot reproduce a nonlinear local grid shift. |
| An Earth transmitter, moving lunar reflector, Moon-orbiting receiver, and camera on a Mars spacecraft | Ephemerides, Earth orientation and polar motion, lunar libration, spacecraft attitude/lever arms, camera calibration, clock/time-scale models | `FrameGraph`, `FramedWorldline`, forward/backward retarded paths and `BearingObservation`. Query each provider at its own transmission/scattering/reception event. Geographic conversion to Moon axes represents the same point; it does not place an Earth station on lunar terrain. |
| Radar through a Martian ice sheet from a descending spacecraft | Body shape and orientation, DEM, ice thickness, frequency-dependent complex permittivity/conductivity, boundary normals and roughness, platform path and clocks | Body-local tangent poses, material propagation, worldlines and custom `ScatterResponse`. Homogeneous-material loss is implemented; layered refraction, anisotropic ice, curved rays and boundary scattering require a supplied propagation/response model. |
| Hundreds of tumbling asteroids reflecting one coherent illuminator into a sparse array | Individual spin/worldpaths and shapes, illumination waveform, phase-linked clocks, polarization and complex surface scattering, visibility, reference delay and coherence windows | `received_bistatic_path`, `coherent_echo`, `integrate_coherent`, channel likelihood and trajectory refinement. All echoes are evaluated at a common receiver event. An unresolved ensemble field can reinforce **or cancel**. Shape points can be framed worldlines; an oblate geodetic API is not a triaxial asteroid shape model. |
| A reflector orbiting one star, another stellar-system receiver, and a nearby millimetre baseline | A common ephemeris/time convention, nearby numerical origins, local differential delay solver, precision budget, physical visibility and propagation model | Hierarchical relative frame poses preserve local baselines beneath huge shared translations. `LorentzBoost` transforms inertial events and velocities when observer simultaneity changes. General relativity, precision apparent astrometry and interstellar medium are external models, not rigid-frame transforms. |
| Arctic RF/polarimetric array correlated with GNSS, ionosondes, balloon profiles and aircraft images during a magnetic storm | Timestamped measurements and uncertainties; GNSS carrier ambiguities/biases; ionosonde response and frequency coverage; balloon trajectories; weather profiles; calibrated images; time-varying magnetic/plasma field and error lineage | `MagneticField`, `ElectronDensity`, frame adapters, field-aligned bases, field-line tracing, TEC/Faraday integration, correlated Gaussian blocks, bearing likelihood and atlas scalar fields. A dipole is a synthetic baseline, not a measured storm model. HF ionosonde rays need magnetoionic ray tracing; high-frequency first-order TEC formulas cannot substitute for it. |
| An incoming X-ray impulse coincident with an optical flash and an RF transient | Exposure/live-time, photon energies/background, detector response, waveform envelopes, arrival-time calibrations, source position/path and association hypotheses | `PhotonObservation`, `BearingObservation`, correlated observables and `JointLikelihood`. Estimate a common source/event through modality-specific likelihoods. Photon counts and unrelated RF carrier phase cannot be summed as one complex field. A burst association must be tested against coincidence/background alternatives. |
| Historical reprocessing of a year of distributed sensing after new clocks, weather or ephemerides become available | Original samples or sufficient statistics, unique observation IDs, original calibration versions, source lineage, time coverage, new nuisance model and process dynamics | Recompute forward observables, profile channel reflectivity, refine local trajectories, update/smooth states and build C1 worldpaths. Six-state recursive tracking is local constant-velocity dynamics; long gravitational arcs and jointly estimated environmental fields require a larger dynamics/hypothesis model. Missing samples/phase/bandwidth cannot be reconstructed from continuity alone. |

## Coordinate conversion order

Normalize units at ingestion. Angle types distinguish radians, degrees, arcseconds and hour angle; feet distinguish international and US survey definitions. Pressure, temperature, frequency, duration and magnetic flux density have explicit conversions. These wrappers validate finite values; signed lengths/durations remain useful, while physically positive frequency/density checks belong to the consuming model. Covariance scales by `J C Jᵀ`, not the factor applied to a mean: metre/foot variance needs the squared length factor.

For a DEM height `H` and source geoid `N`, first form `h = H + N`. `GeographicCoordinate` explicitly follows its `GeodeticReference` height convention; `Geodetic` itself always carries ellipsoidal height. Each geoid must belong to the source/target body's exact horizontal reference. No datum-name registry infers that relationship.

For a rigid horizontal reference change, use `source.convert(point, target, graph, epoch)`. For a scaled/dynamic Helmert change, perform the explicit Cartesian step:

```rust
use geoid_atlas::coordinates::Ellipsoid;
use geoid_atlas::formats::Helmert;
use geoid_atlas::frames::{FrameGraph, State};
use geoid_atlas::reference::*;

let graph = FrameGraph::new("coordinate container");
// Separate source/target descriptions; their actual relation is the supplied Helmert.
let source = GeodeticReference {
    frame: graph.root(), ellipsoid: Ellipsoid::WGS84,
    height: HeightReference::Ellipsoidal,
};
let target = GeodeticReference {
    frame: graph.root(), ellipsoid: Ellipsoid::WGS84,
    height: HeightReference::Ellipsoidal,
};
let parameters = Helmert::from_proj(
    "+proj=helmert +x=0.02 +dx=0.001 +t_epoch=2010"
)?;
let point = GeographicCoordinate {
    latitude_deg: 80.0, longitude_deg: -40.0, height_m: 12.0,
};
let shifted = parameters.at_decimal_year(2026.0)?
    .apply(State::stationary(source.to_cartesian(point)?))?;
let converted = target.from_cartesian(shifted.position_m)?;
assert!(converted.height_m.is_finite());
# Ok::<(), geoid_atlas::Error>(())
```

Those example coefficients are synthetic. A genuine operation uses published coefficients, their validity/uncertainty, and the documented direction/convention. `CartesianTransform::unapply` uses the exact inverse of its linearized Helmert matrix; negating seven coefficients is not that inverse. Its full six-state Jacobian includes scale/rotation rate coupling. Transform uncertainty is conditional on exact parameters; coefficient uncertainty must be propagated separately.

For plate motion, a transform's parameter reference epoch and the coordinate observation epoch have different roles. Transform rates account for reference-frame evolution; they do not automatically propagate the actual station along its tectonic worldline. Propagate the station to the desired event before applying the transformation.

`Planetocentric` is east-positive longitude, centre-ray latitude and radius. `Geodetic` is ellipsoid-normal latitude and height. For west-positive longitude feeds, explicitly negate longitude before ingestion. UTM/Web Mercator implementations are WGS84-specific; arbitrary projected CRSs need a supplied projection/PROJ adapter. Deep body-centre geodetic inverses are deliberately rejected because the nearest ellipsoid-normal coordinates can be ambiguous.

## File-format foothold

`formats::Helmert::from_proj` reads a single Cartesian `+proj=helmert` definition with metre translations, arcsecond rotations, ppm scale, and optional rates per year and `+t_epoch`. Rotations require an explicit `position_vector` or `coordinate_frame` convention. Coordinate input/output are SI Cartesian states. The small-angle matrix matches the non-`+exact` operation; rotations exceeding 1e-3 rad fail. It rejects unknown/duplicate options, nonfinite values, missing rate epochs, `+exact`, `+inv`, pipelines, grids and implicit axis/unit changes. It does not resolve EPSG identities, validate coefficient provenance or parse an entire CRS.

Decimal-year input is a separate convention from `Epoch`: it must match the published coefficient convention. Velocity rate conversion uses a Julian year (31,557,600 s). A civil decimal-year calendar, UTC leaps and barycentric time conversions require explicit adapters.

`read_world_file` reads the six ESRI world-file lines in **A,D,B,E,C,F** order, locating pixel centres, including image rotation/shear. It supplies no CRS, units, height reference, pixels, intrinsics or aircraft pose. Orthorectified DRGs and raw camera pixels need different models. ESRI ASCII DEM ingestion is also implemented. Formats worth adding next are NTv2/NADCON displacement grids, GeoTIFF/COG georeferencing, WKT2/PROJJSON metadata, SPICE ephemeris/orientation adapters, and IONEX/RINEX/SINEX/IERS ingestion with strict epoch, reference and coverage handling. These formats are not currently decoded.

## A polar fusion solve

1. Establish a common coordinate time, inertial frame and calibrated clocks. Enter stations and aircraft/balloons through body/attitude/local-frame adapters. Preserve each measurement's units, uncertainty, original ID and provenance. Do not infer a day-long coverage claim from endpoint timestamps.
2. Build local environmental hypotheses from balloon profiles, weather observations, ionosonde inversions, GNSS slant TEC and magnetometers. GNSS estimates are path integrals; a vertical TEC map does not uniquely specify a 3D electron density. Ionosonde inversion is model-dependent. Weather sea-level pressure is not station pressure.
3. Predict each observation with its response model. RF uses retarded paths, dispersive/group versus carrier phase delay, polarization rotation, clocks and reflectivity. Images use calibrated bearings, emission light time, attitude and exposure. X-rays use photon responses/background and event-time association. The same environmental hypothesis feeds the relevant models.
4. Keep shared nuisances together. `GaussianObservation` accepts an arbitrary full covariance, including mixed physical units. `ObservableModel` supplies the predictions. Put GNSS/weather/magnetometer errors that share calibration or an inversion into one correlated block; independent blocks can enter `JointLikelihood`. Singular covariance fails instead of implying infinite information. The fixed-covariance Gaussian cost omits constants; a varying covariance requires a custom likelihood including its log determinant.
5. Refine only the region/time window relevant to the candidate, checking path/field integration convergence and revisiting uncertainty. Existing trajectory optimization operates on worldlines; simultaneous inference of environmental volumes, clocks and scattering parameters requires a larger optimizer/hypothesis around these forward providers. Conditioning on an environment is not the same as estimating it jointly.

`polar_fusion` implements a synthetic calibration-payload demonstration with rotating Earth, an 80°N station, an aircraft bearing, correlated weather/GNSS observations, a coherent radar reflection, a separate photon counter, a plasma column and field-aligned axes. It evaluates the same physical path across modalities and asserts that a displaced hypothesis costs more. All calibration/ephemeris/environment parameters are fixed, and the photon source is a calibrated test payload. It demonstrates composition, not a measured forecast, storm reconstruction or operational acquisition system.

Magnetic vectors are rotated, not translated like positions. `FramedMagneticField` queries the provider in local coordinates at the event epoch and rotates its result into root axes. The `field_aligned_basis` z-axis follows B, with an explicit transverse guide; a zero field or parallel guide has no unique basis. Geographic longitude is arbitrary at either pole, while field-aligned coordinates remain usable where B is nonzero.

`trace_field_line` integrates `dx/ds = ± B/|B|` with RK4 at a frozen epoch. It can map conjugate hemispheres only when the supplied field and integration coverage support that connection. It is not particle drift, auroral timing, plasma transport or proof of instantaneous north/south correlation. Storm-time topology, reconnection and retarded temporal responses need additional physics.

`integrate_plasma` samples electron density and B along a supplied straight segment at interpolated coordinate times. It returns `∫ne ds` and signed `∫ne B·ds`; the latter gives first-order Faraday rotation at sufficiently high frequency. Reversing the geometric ray reverses the magnetic integral but preserves TEC. For curved/strongly dispersive paths, supply actual ray legs and a coupled propagation solver. Quadrature error, field covariance and ray uncertainty are not included automatically.

## What coherent asteroid recovery actually requires

For receiver time `t_r`, solve each asteroid's scattering event backward, then solve the illuminator's emission event backward. `coherent_echo` applies carrier phase relative to a nearby calibrated emission reference and sums complex fields in **one phase-linked channel**. Its response provider can depend on the actual target spin, bistatic angle, polarization, range losses, medium phase, occultation and waveform envelope. The tests demonstrate quarter-wavelength radial separation causing monostatic cancellation and half-wavelength separation causing reinforcement.

Known trajectories alone do not reveal unknown asteroid surface phase. Model each complex response, profile it in sufficiently short windows, or integrate independent target/channel likelihoods. Natural unrelated asteroids need not remain mutually coherent. Range/Doppler gating, broadband envelopes, receiver arrays and spin diversity can separate contributions; a single ambiguous CW sample may not identify them. Phase focusing should improve a validated estimator against an unfocused baseline rather than presuming every added reflector adds positive power.

Split epochs avoid multiplying a carrier by an enormous absolute timestamp, but do not make arbitrary differential delay exact. At 10 GHz, a 0.1 rad phase budget corresponds to about **1.59 ps** and **0.477 mm** of one-way path. Ephemerides, clocks, atmosphere, target shape and numerical error share that budget. Distant optical interferometry generally needs specialized high-precision differential propagation and source coherence models beyond these vacuum double-precision solvers.

## Relativistic and temporal boundaries

`FrameGraph` performs Newtonian rigid kinematics at one declared coordinate time. `LorentzBoost` separately transforms local inertial events/velocities, preserves the Minkowski interval and changes simultaneity. An array snapshot simultaneous in the source observer generally is not simultaneous in the boosted observer; recompute the worldlines on the desired target-time slice rather than just rotating all its positions. A boost's supplied source-origin event maps to target zero, making both local origins explicit.

Multiple-star frame trees do not themselves supply general relativity. Strong gravity, time-dependent metrics, barycentric clock conventions, gravitational lensing, frequency transport and precision stellar apparent positions require appropriate external models. Sensor proper-time clocks must be connected to coordinate time before using event solvers. The existing point-mass/Shapiro terms are approximations with stated limits, not a complete relativistic propagation engine.

A relativistic boost also mixes electromagnetic fields. `LorentzBoost::electromagnetic_field` transforms E and B together; rotating a magnetic vector alone is appropriate for the rigid frame adapter, not a relativistic observer change. Its tests preserve both electromagnetic invariants.

No adapter silently substitutes a model for absent observations. Providers can return `NoData` or `OutsideCoverage`, and likelihood evaluation propagates those failures. The goal is a composable numerical library whose assumptions can be tested, rather than a label-based promise to translate every possible physical coordinate system without its defining inputs.
