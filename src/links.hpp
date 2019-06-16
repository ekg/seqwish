#ifndef LINKS_HPP_INCLUDED
#define LINKS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include "seqindex.hpp"
#include "mmmultimap.hpp"
#include "pos.hpp"

namespace seqwish {


void derive_links(seqindex_t& seqidx,
                  size_t graph_length,
                  mmmulti::map<uint64_t, pos_t>& path_mm,
                  mmmulti::map<pos_t, pos_t>& link_fwd_mm,
                  mmmulti::map<pos_t, pos_t>& link_rev_mm);

}

#endif
