extern crate cbindgen;

use std::env;
use std::process::Command;

fn main() {
    // Bake `git describe` into the binary for `--version` (never the crate version).
    // Falls through to version.rs's "unknown" if git is unavailable (e.g. crates.io).
    if let Some(git_version) = git_describe() {
        println!("cargo:rustc-env=SEQWISH_GIT_VERSION={git_version}");
        println!("cargo:rerun-if-changed=.git/HEAD");
    }

    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_language(cbindgen::Language::C)
        .with_cpp_compat(true) // Add extern "C" guards for C++
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file("seqwish_rs.h");
}

fn git_describe() -> Option<String> {
    let output = Command::new("git")
        .args(["describe", "--tags", "--long", "--always"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let version = String::from_utf8(output.stdout).ok()?.trim().to_string();
    (!version.is_empty()).then_some(version)
}
