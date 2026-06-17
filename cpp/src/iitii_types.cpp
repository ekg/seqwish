#include "iitii_types.hpp"

namespace seqwish {

uint64_t range_get_beg(const range_pos_t& p) { return p.start; }
uint64_t range_get_end(const range_pos_t& p) { return p.end; }

}
