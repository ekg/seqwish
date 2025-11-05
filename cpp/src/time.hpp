#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>
#include "seqwish_rs.h"

namespace seqwish {

// Thin wrapper around Rust implementation
inline uint64_t time_since_epoch_ms(void) {
    return ::time_since_epoch_ms();
}

// Keep this as inline C++ since it depends on C++ chrono types
inline double seconds_since(const std::chrono::time_point<std::chrono::steady_clock>& then) {
    return ((std::chrono::duration<double>)(std::chrono::steady_clock::now() - then)).count();
}

}
