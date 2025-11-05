// version.hpp: Version reflection information for seqwish builds.
// Modified to use Rust implementation via FFI

#pragma once

#include <string>
#include "seqwish_rs.h"

namespace seqwish {

using namespace std;

/// Class for holding seqwish version - now backed by Rust implementation
class Version {
public:
    /// The Git version description of this build of seqwish
    const static string VERSION;

    /// Get only the version (like v1.7.0-68-g224e7625).
    static string get_version() {
        char* version_c = ::version_get_version();
        if (version_c == nullptr) return "";
        string result(version_c);
        ::version_free_string(version_c);
        return result;
    }

    /// Get the release Git tag version of seqwish that the current version
    /// is based on (e.g. v1.7.0-68-g224e7625 will report v1.7.0).
    static string get_release() {
        char* release_c = ::version_get_release();
        if (release_c == nullptr) return "";
        string result(release_c);
        ::version_free_string(release_c);
        return result;
    }

    /// Get the codename of our released version
    static string get_codename() {
        char* codename_c = ::version_get_codename();
        if (codename_c == nullptr) return "";
        string result(codename_c);
        ::version_free_string(codename_c);
        return result;
    }

    /// Get a short one-line description of the current version with no terminating newline.
    static string get_short() {
        char* short_c = ::version_get_short();
        if (short_c == nullptr) return "";
        string result(short_c);
        ::version_free_string(short_c);
        return result;
    }

private:
    // Not constructable
    Version() = delete;
};

}
