#pragma once

#include "mmmultimap.hpp"
#include "mmiitree.hpp"
#include "pos.hpp"

namespace seqwish {

struct match_t {
    uint64_t start = 0; // from q start
    uint64_t end = 0;   // to q end
    pos_t pos = 0;      // where it matches
};

match_t get_match(mmmulti::iitree<uint64_t, pos_t>& iitree, uint64_t idx);

bool operator<(const match_t& a, const match_t& b);

bool operator==(const match_t& a, const match_t& b);

}
