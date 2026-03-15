use std::env;
use std::path::PathBuf;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    let config = cbindgen::Config::from_file(format!("{}/cbindgen.toml", crate_dir))
        .expect("failed to read cbindgen.toml");

    let header = cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
        .expect("failed to generate C bindings");

    // Write to the cargo output directory.
    let out_header = out_dir.join("mordred.h");
    header.write_to_file(&out_header);

    // Also copy to mordred-ffi/include/ so downstream consumers can find it
    // without digging into target/.
    let include_dir = PathBuf::from(&crate_dir).join("include");
    std::fs::create_dir_all(&include_dir).ok();
    let include_header = include_dir.join("mordred.h");
    header.write_to_file(&include_header);
}
