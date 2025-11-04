#include "version.hpp"

// Functions now implemented in Rust (seqwish-rs/src/version.rs)
// This file only provides the static VERSION string initialization

namespace seqwish {

// Initialize the static VERSION string by calling the Rust implementation
const string Version::VERSION = Version::get_version();

}
