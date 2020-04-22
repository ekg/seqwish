#include "time.hpp"

namespace seqwish {

uint64_t time_since_epoch_ms(void) {
  using namespace std::chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

double seconds_since(const std::chrono::time_point<std::chrono::steady_clock>& then) {
    return ((std::chrono::duration<double>)(std::chrono::steady_clock::now() - then)).count();
}

}
