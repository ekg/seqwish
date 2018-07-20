#include "pos.hpp"

namespace seqwish {

bool operator<(const aln_pos_t& a, const aln_pos_t& b) {
    return a.pos < b.pos && a.aln_length < b.aln_length;
}

bool operator==(const aln_pos_t& a, const aln_pos_t& b) {
    return a.pos == b.pos && a.aln_length == b.aln_length;
}

pos_t make_pos_t(uint64_t offset, bool is_rev) {
    // top bit is reserved for is_rev flag
    // the rest is our offset in the input sequence vector
    uint64_t rev_mask = (uint64_t)1; // the bit mask
    pos_t pos = offset<<1;
    // https://graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    pos = (pos & ~rev_mask) | (-is_rev & rev_mask);
    return pos;
}

uint64_t offset(const pos_t& pos) {
    //return (pos & ~(uint64_t)1) >> 1;
    return pos >> 1;
}

bool is_rev(const pos_t& pos) {
    return pos & (uint64_t)1;
}

void incr_pos(pos_t& pos) {
    if (is_rev(pos)) {
        pos -= 2;
    } else {
        pos += 2;
    }
}

void incr_pos(pos_t& pos, size_t by) {
    if (is_rev(pos)) {
        pos -= 2*by;
    } else {
        pos += 2*by;
    }
}

pos_t rev_pos_t(const pos_t& pos) {
    return make_pos_t(offset(pos), !is_rev(pos));
}

std::string pos_to_string(const pos_t& pos) {
    return std::to_string(offset(pos)) + (is_rev(pos)?"-":"+");
}

}
