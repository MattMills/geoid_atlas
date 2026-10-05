//! Ground-path primitives for an FT8/interference port. Heights and propagation
//! conditions are separate inputs; locator centres below are representative sites.
use geoid_atlas::Result;
use geoid_atlas::coordinates::Geodetic;
use geoid_atlas::geodesics::{SurfaceGeodesic, great_circle_m, ground_distance_m};
use geoid_atlas::maidenhead::MaidenheadCell;

fn main() -> Result<()> {
    let cell: MaidenheadCell = "IO91wm".parse()?;
    let start = cell.centre();
    println!("{} cell bounds: {:?}", cell.locator(), cell.bounds());
    println!(
        "Uniform-cell tangent uncertainty: {:?}",
        cell.uniform_uncertainty()
    );
    for (name, end) in [
        ("Sydney", Geodetic::new(-33.9, 151.2, 0.0)?),
        ("Wellington", Geodetic::new(-41.3, 174.78, 0.0)?),
    ] {
        let path = SurfaceGeodesic::wgs84().path(start, end);
        let exact = ground_distance_m(start, end);
        let sphere = great_circle_m(start, end);
        let midpoint = path.midpoint();
        assert!((ground_distance_m(start, midpoint) - exact / 2.0).abs() < 1e-7);
        println!(
            "{name}: WGS84={exact:.3} m, mean-radius sphere={sphere:.3} m, difference={:.3} m",
            exact - sphere
        );
        println!(
            "  forward headings: {:.6} -> {:.6} deg; surface midpoint: {:?}",
            path.inverse().initial_azimuth_deg,
            path.inverse().final_azimuth_deg,
            midpoint
        );
        let radii = midpoint.curvature_radii();
        let heading = path.position(exact / 2.0)?.final_azimuth_deg;
        println!(
            "  local midpoint curvature: M={:.3} m, N={:.3} m, along path={:.3} m",
            radii.meridional_m(),
            radii.prime_vertical_m(),
            radii.along_azimuth_m(heading)?
        );
    }
    Ok(())
}
