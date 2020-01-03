#include "match.hpp"

namespace seqwish {

match_t get_match(mmmulti::iitree<uint64_t, pos_t>& iitree, uint64_t idx) {
    return { iitree.start(idx), iitree.end(idx), iitree.data(idx) };
}

bool operator<(const match_t& a, const match_t& b) {
    return a.start < b.start && a.end < b.end && a.pos < b.pos;
}

bool operator==(const match_t& a, const match_t& b) {
    return a.start == b.start && a.end == b.end && a.pos == b.pos;
}

}
