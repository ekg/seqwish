#pragma once

#include "iitii.h"
#include "pos.hpp"

namespace seqwish {

struct range_pos_t {
    uint64_t start;
    uint64_t end;
    pos_t pos;
};

inline uint64_t range_get_beg(const range_pos_t& p) { return p.start; }
inline uint64_t range_get_end(const range_pos_t& p) { return p.end; }
using range_pos_iitii = iitii::iitii<uint64_t, range_pos_t, range_get_beg, range_get_end>;

}
