#include "pos.hpp"

namespace seqwish {

pos_t make_pos_t(uint64_t offset, bool is_rev) {
    // top bit is reserved for is_rev flag
    // the rest is our offset in the input sequence vector
    uint64_t rev_mask = (uint64_t)1<<63; // the bit mask
    pos_t pos;
    // https://graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    pos = (pos & ~rev_mask) | (-is_rev & rev_mask);
    return pos;
}

uint64_t offset(const pos_t& pos) {
    return pos & ~(uint64_t)1<<63;
}

bool is_rev(const pos_t& pos) {
    return pos & (uint64_t)1<<63;
}

}
