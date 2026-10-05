fn main() {
    println!("cargo:rerun-if-env-changed=GEOID_NATIVE_LIB_DIR");
    if let Some(directory) = std::env::var_os("GEOID_NATIVE_LIB_DIR") {
        println!(
            "cargo:rustc-link-search=native={}",
            std::path::Path::new(&directory).display()
        );
    }
    for (feature, library) in [
        ("CARGO_FEATURE_PROJ_BACKEND", "proj"),
        ("CARGO_FEATURE_GDAL_BACKEND", "gdal"),
    ] {
        if std::env::var_os(feature).is_some() {
            println!("cargo:rustc-link-lib={library}");
        }
    }
}
