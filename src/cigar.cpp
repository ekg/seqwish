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

}
