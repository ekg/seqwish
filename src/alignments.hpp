#ifndef ALIGNMENTS_HPP_INCLUDED
#define ALIGNMENTS_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include "paf.hpp"
#include "sxs.hpp"
#include "mmmultimap.hpp"
#include "seqindex.hpp"
#include "gzstream.h"
#include "pos.hpp"

namespace seqwish {


void unpack_paf_alignments(const std::string& paf_file,
                           mmmulti::map<pos_t, pos_t>& aln_mm,
                           seqindex_t& seqidx);

void unpack_sxs_alignments(const std::string& sxs_file,
                           mmmulti::map<pos_t, pos_t>& aln_mm,
                           seqindex_t& seqidx);

/*
void filter_alignments(mmmulti::map<pos_t, aln_pos_t>& aln_mm,
                       mmmulti::map<pos_t, pos_t>& aln_filt_mm,
                       uint64_t aln_min_length,
                       uint64_t aln_keep_n_longest,
                       seqindex_t& seqidx);
*/

}

#endif
