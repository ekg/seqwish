#pragma once

#include <cstdint>
#include <string>
#include "seqwish_rs.h"

namespace seqwish {

typedef uint64_t pos_t;

struct aln_pos_t {
    pos_t pos;
    uint64_t aln_length;
};

// Operators for aln_pos_t
inline bool operator<(const aln_pos_t& a, const aln_pos_t& b) {
    return a.pos < b.pos && a.aln_length < b.aln_length;
}

inline bool operator==(const aln_pos_t& a, const aln_pos_t& b) {
    return a.pos == b.pos && a.aln_length == b.aln_length;
}

// Position manipulation functions - now thin wrappers around Rust implementation
inline pos_t make_pos_t(uint64_t offset, bool is_rev) {
    return pos_make_pos_t(offset, is_rev);
}

inline uint64_t offset(const pos_t& pos) {
    return pos_offset(pos);
}

inline bool is_rev(const pos_t& pos) {
    return pos_is_rev(pos);
}

inline void incr_pos(pos_t& pos) {
    pos_incr_pos(&pos);
}

inline void incr_pos(pos_t& pos, size_t by) {
    pos_incr_pos_by(&pos, by);
}

inline void decr_pos(pos_t& pos) {
    pos_decr_pos(&pos);
}

inline void decr_pos(pos_t& pos, size_t by) {
    pos_decr_pos_by(&pos, by);
}

inline pos_t rev_pos_t(const pos_t& pos) {
    return pos_rev_pos_t(pos);
}

inline std::string pos_to_string(const pos_t& pos) {
    char* result = pos_to_string_c(pos);
    if (result == nullptr) {
        return "";
    }
    std::string str(result);
    temp_file_free_string(result);  // Reuse the free function
    return str;
}

}
