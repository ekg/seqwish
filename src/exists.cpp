#include "exists.hpp"
#include "seqwish_rs.h"

namespace seqwish {

bool file_exists(const std::string& name) {
    return ::file_exists(name.c_str());
}

}
