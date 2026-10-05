use geoid_atlas::coordinates::{Ellipsoid, Geodetic};
use geoid_atlas::maidenhead::MaidenheadCell;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn standard_locator_centres_and_cell_dimensions() {
    for (locator, lon, lat, lon_width, lat_width) in [
        ("JO", 10.0, 55.0, 20.0, 10.0),
        ("JO22", 5.0, 52.5, 2.0, 1.0),
        ("IO91wm", -0.125, 51.52083333333333, 1.0 / 12.0, 1.0 / 24.0),
        (
            "FN31pr",
            -72.70833333333333,
            41.72916666666667,
            1.0 / 12.0,
            1.0 / 24.0,
        ),
        ("FN31", -73.0, 41.5, 2.0, 1.0),
        (
            "AA00aa00",
            -180.0 + 1.0 / 240.0,
            -90.0 + 1.0 / 480.0,
            1.0 / 120.0,
            1.0 / 240.0,
        ),
    ] {
        let cell: MaidenheadCell = locator.parse().unwrap();
        let centre = cell.centre();
        close(centre.longitude_deg(), lon, 1e-12);
        close(centre.latitude_deg(), lat, 1e-12);
        close(cell.bounds().longitude_width_deg(), lon_width, 1e-12);
        close(cell.bounds().latitude_width_deg(), lat_width, 1e-12);
        assert_eq!(cell.locator(), locator.to_ascii_uppercase());
        assert_eq!(centre.height_m(), 0.0);
        assert_eq!(cell.centre_at_height(-150.0).unwrap().height_m(), -150.0);
        assert!(cell.centre_at_height(f64::NAN).is_err());
    }
}

#[test]
fn precision_levels_and_world_boundaries() {
    let mut longitude_width = 20.0;
    let mut latitude_width = 10.0;
    for (i, locator) in [
        "RR",
        "RR99",
        "RR99XX",
        "RR99XX99",
        "RR99XX99XX",
        "RR99XX99XX99",
    ]
    .iter()
    .enumerate()
    {
        if i > 0 {
            let radix = if i % 2 == 1 { 10.0 } else { 24.0 };
            longitude_width /= radix;
            latitude_width /= radix;
        }
        let cell = MaidenheadCell::parse(locator).unwrap();
        let bounds = cell.bounds();
        close(bounds.east_deg, 180.0, 1e-12);
        close(bounds.north_deg, 90.0, 1e-12);
        close(bounds.west_deg, 180.0 - longitude_width, 1e-12);
        close(bounds.south_deg, 90.0 - latitude_width, 1e-12);
        let centre = cell.centre();
        assert!(centre.latitude_deg() < 90.0 && centre.longitude_deg() < 180.0);
    }
    let origin = MaidenheadCell::parse("AA00AA00AA00").unwrap().bounds();
    assert_eq!(origin.west_deg, -180.0);
    assert_eq!(origin.south_deg, -90.0);
    // Neighbouring cells meet exactly to numerical precision without overlap.
    let west = MaidenheadCell::parse("IO91WM").unwrap().bounds();
    let east = MaidenheadCell::parse("IO91XM").unwrap().bounds();
    close(west.east_deg, east.west_deg, 1e-12);
    let north = MaidenheadCell::parse("IO91WN").unwrap().bounds();
    close(west.north_deg, north.south_deg, 1e-12);
}

#[test]
fn cell_uncertainty_states_uniform_assumption_and_local_scale() {
    let cell = MaidenheadCell::parse("JJ00AA").unwrap();
    let centre = cell.centre();
    let bounds = cell.bounds();
    let uncertainty = cell.uniform_uncertainty();
    close(
        uncertainty.longitude_std_deg,
        bounds.longitude_width_deg() / 12.0_f64.sqrt(),
        1e-15,
    );
    close(
        uncertainty.latitude_std_deg,
        bounds.latitude_width_deg() / 12.0_f64.sqrt(),
        1e-15,
    );
    close(
        uncertainty.east_std_m,
        uncertainty.east_width_m / 12.0_f64.sqrt(),
        1e-9,
    );
    close(
        uncertainty.north_std_m,
        uncertainty.north_width_m / 12.0_f64.sqrt(),
        1e-9,
    );
    // Independent Cartesian finite differences approximate the cell's tangent size.
    let east = Geodetic::new(centre.latitude_deg(), bounds.east_deg, 0.0)
        .unwrap()
        .to_ecef();
    let west = Geodetic::new(centre.latitude_deg(), bounds.west_deg, 0.0)
        .unwrap()
        .to_ecef();
    let north = Geodetic::new(bounds.north_deg, centre.longitude_deg(), 0.0)
        .unwrap()
        .to_ecef();
    let south = Geodetic::new(bounds.south_deg, centre.longitude_deg(), 0.0)
        .unwrap()
        .to_ecef();
    close(uncertainty.east_width_m, east.distance(west), 0.001);
    close(uncertainty.north_width_m, north.distance(south), 0.001);
    let small = cell
        .uniform_uncertainty_on(Ellipsoid::new(1_000.0, 0.0).unwrap())
        .unwrap();
    close(
        small.north_width_m,
        1_000.0 * bounds.latitude_width_deg().to_radians(),
        1e-12,
    );
    close(
        small.east_width_m,
        1_000.0
            * centre.latitude_deg().to_radians().cos()
            * bounds.longitude_width_deg().to_radians(),
        1e-12,
    );
    let polar = MaidenheadCell::parse("RR99XX99XX99")
        .unwrap()
        .uniform_uncertainty();
    assert!(polar.east_width_m > 0.0 && polar.north_width_m > 0.0);
    assert!(polar.east_width_m < polar.north_width_m);
}

#[test]
fn strict_locator_validation() {
    for invalid in [
        "",
        "A",
        "AAA",
        "AA0",
        "AA00A",
        "AA00AA0",
        "AA00AA00AA0",
        "SA",
        "AZ",
        "00",
        "AAZZ",
        "AA0A",
        "AA00YY",
        "AA00AZ",
        "AA00AA0X",
        "AA00AA00YY",
        "AA00AA00AA0A",
        "AA00AA00AA00AA",
        " IO91WM",
        "IO91WM ",
        "IÓ91WM",
        "IO91-WM",
        "AA00\0A",
    ] {
        assert!(
            MaidenheadCell::parse(invalid).is_err(),
            "accepted {invalid:?}"
        );
    }
    // Every field, including the polar and date-line cells, is decodable.
    for lon in b'A'..=b'R' {
        for lat in b'A'..=b'R' {
            let locator = String::from_utf8(vec![lon, lat]).unwrap();
            let centre = MaidenheadCell::parse(&locator).unwrap().centre();
            assert!((-90.0..90.0).contains(&centre.latitude_deg()));
            assert!((-180.0..180.0).contains(&centre.longitude_deg()));
        }
    }
}
