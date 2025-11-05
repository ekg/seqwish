#pragma once

#include <string>
#include <cstdio>
#include "seqwish_rs.h"

/**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 *
 * NOTE: This is now a thin wrapper around the Rust implementation.
 */
namespace temp_file {

    /// Create a temporary file starting with the given base name
    inline std::string create(const std::string& base, const std::string& suffix) {
        char* result = temp_file_create(base.c_str(), suffix.c_str());
        if (result == nullptr) {
            return "";
        }
        std::string path(result);
        temp_file_free_string(result);
        return path;
    }

    /// Remove a temporary file
    inline void remove(const std::string& filename) {
        temp_file_remove(filename.c_str());
    }

    /// Set a temp dir, overriding system defaults and environment variables.
    inline void set_dir(const std::string& new_temp_dir) {
        temp_file_set_dir(new_temp_dir.c_str());
    }

    /// Get the current temp dir
    inline std::string get_dir() {
        char* result = temp_file_get_dir();
        if (result == nullptr) {
            return "";
        }
        std::string dir(result);
        temp_file_free_string(result);
        return dir;
    }

    inline void set_keep_temp(bool setting) {
        temp_file_set_keep_temp(setting);
    }

} // namespace temp_file