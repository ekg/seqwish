#ifndef LINKS_HPP_INCLUDED
#define LINKS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include "seqindex.hpp"
#include "dmultimap.hpp"
#include "pos.hpp"

namespace seqwish {

void derive_links(seqindex_t& seqidx, size_t graph_length, dmultimap<uint64_t, pos_t>& path_mm, dmultimap<pos_t, pos_t>& link_mm);

}

#endif
