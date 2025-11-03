#ifndef CIGAR_HPP_INCLUDED
#define CIGAR_HPP_INCLUDED

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdint>
#include "seqwish_rs.h"

namespace seqwish {

struct cigar_op_t { uint64_t len; char op; };
typedef std::vector<cigar_op_t> cigar_t;

inline cigar_t cigar_from_string(const std::string& s) {
    CigarHandle* handle = ::cigar_from_string(s.c_str());
    if (handle == nullptr) {
        return cigar_t();
    }

    cigar_t result;
    size_t len = ::cigar_length(handle);
    result.reserve(len);

    for (size_t i = 0; i < len; ++i) {
        uint64_t op_len;
        uint8_t op;
        if (::cigar_get_op(handle, i, &op_len, &op)) {
            result.push_back({op_len, static_cast<char>(op)});
        }
    }

    ::cigar_free(handle);
    return result;
}

inline std::string cigar_to_string(const cigar_t& cigar) {
    // Build C-compatible vector
    std::vector<uint64_t> lengths;
    std::vector<uint8_t> ops;
    lengths.reserve(cigar.size());
    ops.reserve(cigar.size());

    for (const auto& op : cigar) {
        lengths.push_back(op.len);
        ops.push_back(static_cast<uint8_t>(op.op));
    }

    // Build string directly (simpler than converting through handle)
    std::stringstream ss;
    for (const auto& elem : cigar) {
        ss << elem.len << elem.op;
    }
    return ss.str();
}

}

#endif
