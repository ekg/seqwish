#ifndef ALIGNMENTS_HPP_INCLUDED
#define ALIGNMENTS_HPP_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include "paf.hpp"
#include "sxs.hpp"
#include "mmmultimap.hpp"
#include "mmiitree.hpp"
#include "iitii_types.hpp"
#include "seqindex.hpp"
#include "gzstream.h"
#include "pos.hpp"

namespace seqwish {


void unpack_paf_alignments(const std::string& paf_file,
                           range_pos_iitii::builder& aln_iitree_builder,
                           seqindex_t& seqidx,
                           uint64_t min_match_len);

void unpack_sxs_alignments(const std::string& sxs_file,
                           range_pos_iitii::builder& aln_iitree_builder,
                           seqindex_t& seqidx,
                           uint64_t min_match_len);

/*
void filter_alignments(mmmulti::map<pos_t, aln_pos_t>& aln_mm,
                       mmmulti::map<pos_t, pos_t>& aln_filt_mm,
                       uint64_t aln_min_length,
                       uint64_t aln_keep_n_longest,
                       seqindex_t& seqidx);
*/

}

#endif
