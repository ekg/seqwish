#pragma once

#include <string>
#include <mutex>
#include <cstdio>
#include <dirent.h>

#include <cstring>

/**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 */
namespace temp_file {

    /// Create a temporary file starting with the given base name
    std::string create(const std::string& base, const std::string& suffix);

    /// Remove a temporary file
    void remove(const std::string& filename);

    /// Set a temp dir, overriding system defaults and environment variables.
    void set_dir(const std::string& new_temp_dir);

    /// Get the current temp dir
    std::string get_dir();

    void set_keep_temp(bool setting);

} // namespace temp_file