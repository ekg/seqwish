#ifndef DNA_HPP_INCLUDED
#define DNA_HPP_INCLUDED

#include <string>
#include "seqwish_rs.h"

namespace seqwish {

// Note: Despite the name, this is just complement, not reverse complement
inline char dna_reverse_complement(const char& c) {
    return ::dna_complement(static_cast<uint8_t>(c));
}

inline std::string dna_reverse_complement(const std::string& seq) {
    std::string result(seq.size(), '\0');
    ::dna_reverse_complement(seq.data(), seq.size(), &result[0]);
    return result;
}

inline void dna_reverse_complement_in_place(std::string& seq) {
    ::dna_reverse_complement_in_place(&seq[0], seq.size());
}

}

#endif
