#include "time.hpp"
#include "seqwish_rs.h"

namespace seqwish {

uint64_t time_since_epoch_ms(void) {
  return ::time_since_epoch_ms();
}

double seconds_since(const std::chrono::time_point<std::chrono::steady_clock>& then) {
    return ((std::chrono::duration<double>)(std::chrono::steady_clock::now() - then)).count();
}

}
