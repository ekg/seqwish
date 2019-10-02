#include "cigar.hpp"

namespace seqwish {

cigar_t cigar_from_string(const std::string& s) {
    cigar_t cigar;
    std::string number;
    char type = 0;
    for (auto c : s) {
        if (isdigit(c)) {
            if (type == 0) {
                number += c;
            } else {
                cigar.push_back({std::stoul(number), type});
                number.clear();
                type = 0;
                number += c;
            }
        } else {
            type = c;
        }
    }
    if (!number.empty() && type != 0) {
        cigar.push_back({std::stoul(number), type});
    }
    return cigar;
}

std::string cigar_to_string(const cigar_t& cigar) {
    std::stringstream ss;
    for (auto& elem : cigar) {
        ss << elem.len << elem.op;
    }
    return ss.str();
}

size_t cigar_query_length(const cigar_t& cigar) {
    size_t length = 0;
    for (auto& elem : cigar) {
        switch (elem.op) {
        case 'M':
        case 'I':
            length += elem.len;
            break;
        default:
            break;
        }
    }
    return length;
}

size_t cigar_target_length(const cigar_t& cigar) {
    size_t length = 0;
    for (auto& elem : cigar) {
        switch (elem.op) {
        case 'M':
        case 'D':
            length += elem.len;
            break;
        default:
            break;
        }
    }
    return length;
}

}
