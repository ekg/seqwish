#pragma once

#include <cstdint>
#include <string>

namespace seqwish {

typedef uint64_t pos_t;
struct aln_pos_t { pos_t pos; uint64_t aln_length; };
bool operator<(const aln_pos_t& a, const aln_pos_t& b);
bool operator==(const aln_pos_t& a, const aln_pos_t& b);
pos_t make_pos_t(uint64_t offset, bool is_rev);
uint64_t offset(const pos_t& pos);
bool is_rev(const pos_t& pos);
void incr_pos(pos_t& pos);
void incr_pos(pos_t& pos, size_t by);
void decr_pos(pos_t& pos);
void decr_pos(pos_t& pos, size_t by);
pos_t rev_pos_t(const pos_t& pos);
std::string pos_to_string(const pos_t& pos);

}
