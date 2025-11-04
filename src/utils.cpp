#include <string>
#include "utils.hpp"
#include "seqwish_rs.h"

namespace seqwish {
    double handy_parameter(const std::string& value, const double default_value) {
        return ::handy_parameter(value.c_str(), default_value);
    }
}
