#ifndef EXISTS_HPP_INCLUDED
#define EXISTS_HPP_INCLUDED

#include <string>
#include <sys/stat.h>
#include "seqwish_rs.h"

namespace seqwish {

inline bool file_exists(const std::string& name) {
    return ::file_exists(name.c_str());
}

}

#endif
