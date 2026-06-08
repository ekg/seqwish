#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>

namespace seqwish {

uint64_t time_since_epoch_ms(void);
double seconds_since(const std::chrono::time_point<std::chrono::steady_clock>& then);

}
