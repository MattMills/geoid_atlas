//! Uses the nominal frequency plan in the supplied KiwiSDR solve reports.
//! Optionally audits the rounded geometry table from a local synthesis.txt.
//! Does not fit clocks, delays or cosmological coefficients from summary values.
use geoid_atlas::Result;
use geoid_atlas::observability::*;

fn main() -> Result<()> {
    let f = [
        78_000, 162_000, 540_000, 693_000, 882_000, 954_000, 1_089_000, 1_278_000,
    ];
    let frequencies: Vec<_> = f.iter().map(|f| *f as f64).collect();
    let spectral = clock_dispersion_separability(&frequencies, None)?;
    let resolution = delay_resolution(461.0, DistanceConvention::OneWayPath)?;
    let lattice = nominal_frequency_lattice(&f)?;
    println!(
        "Unweighted f versus 1/f correlation: {:.6}; variance inflation: {:.4}",
        spectral.correlation,
        spectral.variance_inflation.unwrap()
    );
    println!(
        "One-way 461 Hz delay cell: {:.3} km (not a standard deviation)",
        resolution.distance_cell_m / 1000.0
    );
    println!(
        "Nominal pair-beat maximum: {:.3} km; combined relative-phase period: {:.3} km",
        lattice.largest_pair_beat_path_m / 1000.0,
        lattice.one_way_path_period_m / 1000.0
    );
    println!(
        "Lattice values require mutually calibrated carrier phases; the reports do not establish that."
    );
    println!(
        "34 receiver scales/offsets + 8 transmitter terms - 2 gauge freedoms = 74 parameters; 262 cells imply 188 dof if full rank."
    );
    if let Some(path) = std::env::args().nth(1) {
        let text =
            std::fs::read_to_string(path).map_err(|e| geoid_atlas::Error::Parse(e.to_string()))?;
        let section = text
            .split_once("6. THE GEOMETRY AS A TARGET")
            .and_then(|(_, tail)| {
                tail.split_once("7. THE FRAME OVER THE CAPTURE")
                    .map(|(section, _)| section)
            })
            .ok_or_else(|| {
                geoid_atlas::Error::Parse("expected synthesis geometry section".into())
            })?;
        let mut residuals = vec![];
        let (mut inside, mut outside, mut unknown) = (0, 0, 0);
        for line in section.lines() {
            let p: Vec<_> = line.split_whitespace().collect();
            if p.len() < 9 || p[2].parse::<u64>().is_err() {
                continue;
            }
            let (Ok(measured), Ok(geometry), Ok(residual)) = (
                p[4].parse::<f64>(),
                p[5].parse::<f64>(),
                p[6].parse::<f64>(),
            ) else {
                continue;
            };
            if !measured.is_finite() || !geometry.is_finite() || !residual.is_finite() {
                return Err(geoid_atlas::Error::Parse(
                    "nonfinite geometry summary".into(),
                ));
            }
            residuals.push(residual.abs());
            if p[7] == "-" {
                unknown += 1;
            } else {
                let cell = p[7]
                    .parse::<f64>()
                    .map_err(|_| geoid_atlas::Error::Parse("invalid cell width".into()))?;
                if !cell.is_finite() || cell <= 0.0 {
                    return Err(geoid_atlas::Error::Parse("invalid cell width".into()));
                }
                if residual.abs() <= cell {
                    inside += 1;
                } else {
                    outside += 1;
                }
            }
        }
        if residuals.is_empty() {
            return Err(geoid_atlas::Error::NoData);
        }
        residuals.sort_by(f64::total_cmp);
        let n = residuals.len();
        let median = if n % 2 == 0 {
            (residuals[n / 2 - 1] + residuals[n / 2]) / 2.0
        } else {
            residuals[n / 2]
        };
        println!(
            "Local report audit: {n} rounded geometry rows; {inside} inside their reported cell, {outside} outside, {unknown} with no cell; median |reported residual| {median:.1} km"
        );
        println!(
            "These rows lack lag uncertainties, per-slot times and path/clock covariance; no Gaussian refit is performed."
        );
    }
    Ok(())
}
