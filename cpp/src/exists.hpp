#ifndef EXISTS_HPP_INCLUDED
#define EXISTS_HPP_INCLUDED

#include <string>
#include <sys/stat.h>

namespace seqwish {

bool file_exists(const std::string& name);

}

#endif
