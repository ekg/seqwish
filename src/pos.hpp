#ifndef POS_H_INCLUDED
#define POS_H_INCLUDED

#include <cstdint>
#include <string>

namespace seqwish {

typedef uint64_t pos_t;
pos_t make_pos_t(uint64_t offset, bool is_rev);
uint64_t offset(const pos_t& pos);
bool is_rev(const pos_t& pos);
void incr_pos(pos_t& pos);
std::string pos_to_string(const pos_t& pos);

}

#endif
