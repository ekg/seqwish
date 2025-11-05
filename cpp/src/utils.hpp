#pragma once

#include <math.h>
#include <cstdint>
#include <string>
#include "seqwish_rs.h"

namespace seqwish {
    inline double handy_parameter(const std::string& value, const double default_value) {
        return ::handy_parameter(value.c_str(), default_value);
    }
}
