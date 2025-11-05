#include "match.hpp"

namespace seqwish {

match_t get_match(mmmulti::iitree<uint64_t, pos_t>& iitree, uint64_t idx) {
    return { iitree.start(idx), iitree.end(idx), iitree.data(idx) };
}

}
