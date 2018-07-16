#ifndef POS_H_INCLUDED
#define POS_H_INCLUDED

#include <cstdint>

namespace seqwish {

typedef uint64_t pos_t;
pos_t make_pos_t(uint64_t offset, bool is_rev);
uint64_t offset(const pos_t& pos);
bool is_rev(const pos_t& pos);

}

#endif
