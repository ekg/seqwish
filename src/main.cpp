#include <iostream>
#include <cstdint>
#include <cstdio>
#include <string>
#include "args.hxx"
#include "dmultimap.hpp"
#include "sdsl/bit_vectors.hpp"
#include "seqindex.hpp"
#include "paf.hpp"
#include "alignments.hpp"
#include "transclosure.hpp"
#include "links.hpp"
#include "compact.hpp"
#include "gfa.hpp"
#include "vgp.hpp"
#include "pos.hpp"
#include "threads.hpp"
#include "exists.hpp"

using namespace seqwish;

int main(int argc, char** argv) {
    args::ArgumentParser parser("seqwish: a variation graph inducer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> sxs_alns(parser, "FILE", "induce the graph from these .sxs formatted alignments", {'a', "sxs-alns"});
    args::ValueFlag<std::string> paf_alns(parser, "FILE", "induce the graph from these PAF formatted alignments", {'p', "paf-alns"});
    args::ValueFlag<std::string> seqs(parser, "FILE", "the sequences used to generate the alignments (FASTA, FASTQ, .seq)", {'s', "seqs"});
    args::ValueFlag<std::string> base(parser, "BASE", "build graph using this basename", {'b', "base"});
    args::ValueFlag<std::string> gfa_out(parser, "FILE", "write the graph in GFA to FILE", {'g', "gfa"});
    args::ValueFlag<std::string> sml_in(parser, "FILE", "use the sequence match list in FILE to subset the input alignments", {'m', "match-list"});
    args::ValueFlag<std::string> vgp_base(parser, "BASE", "write the graph in VGP format with basename FILE", {'o', "vgp-out"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<uint64_t> repeat_max(parser, "N", "limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during graph induction", {'k', "keep-temp"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    size_t n_threads = args::get(num_threads);
    if (n_threads) {
        omp_set_num_threads(args::get(num_threads));
    } else {
        omp_set_num_threads(1);
    }

    if (!args::get(seqs).empty() && !file_exists(args::get(seqs))) {
        std::cerr << "[seqwish] ERROR: input sequence file " << args::get(seqs) << " does not exist" << std::endl;
        return 2;
    }
    if (!args::get(sxs_alns).empty() && !file_exists(args::get(sxs_alns))) {
        std::cerr << "[seqwish] ERROR: input alignment file " << args::get(sxs_alns) << " does not exist" << std::endl;
        return 4;
    }
    if (!args::get(paf_alns).empty() && !file_exists(args::get(paf_alns))) {
        std::cerr << "[seqwish] ERROR: input alignment file " << args::get(paf_alns) << " does not exist" << std::endl;
        return 4;
    }

    // 1) index the queries (Q) to provide sequence name to position and position to sequence name mapping, generating a CSA and a sequence file
    seqindex_t seqidx;
    seqidx.build_index(args::get(seqs), args::get(base));
    seqidx.save();

    // 2) parse the alignments into position pairs and index (A)
    std::string aln_idx = args::get(base) + ".sqa";
    std::remove(aln_idx.c_str());
    dmultimap<uint64_t, pos_t> aln_mm(aln_idx);
    if (!args::get(sxs_alns).empty()) {
        unpack_sxs_alignments(args::get(sxs_alns), aln_mm, seqidx); // yields array A
        if (args::get(debug)) dump_sxs_alignments(args::get(sxs_alns));
    } else if (!args::get(paf_alns).empty()) {
        unpack_paf_alignments(args::get(paf_alns), aln_mm, seqidx); // yields array A
        if (args::get(debug)) dump_paf_alignments(args::get(paf_alns));
    } else {
        assert(false);
    }
    if (args::get(debug)) {
        aln_mm.for_each_pair([&](const uint64_t& p1, const pos_t& p2) {
                std::cout << "aln_mm" << "\t" << p1 << "\t" << offset(p2) << "\t" << (is_rev(p2)?"-":"+") << std::endl; });
    }

    // 3) find the transitive closures via the alignments and construct S, N, and P indexed arrays
    std::string seq_v_file = args::get(base) + ".sqs";
    std::string node_mm_idx = args::get(base) + ".sqn";
    std::string path_mm_idx = args::get(base) + ".sqp";
    std::remove(seq_v_file.c_str());
    std::remove(node_mm_idx.c_str());
    std::remove(path_mm_idx.c_str());
    dmultimap<uint64_t, pos_t> node_mm(node_mm_idx);
    dmultimap<uint64_t, pos_t> path_mm(path_mm_idx);
    size_t graph_length = compute_transitive_closures(seqidx, aln_mm, seq_v_file, node_mm, path_mm, args::get(repeat_max));
    if (args::get(debug)) {
        node_mm.for_each_pair([&](const uint64_t& p1, const pos_t& p2) {
                std::cout << "node_mm" << "\t" << p1 << "\t" << offset(p2) << "\t" << (is_rev(p2)?"-":"+") << std::endl; });
        path_mm.for_each_pair([&](const uint64_t& p1, const pos_t& p2) {
                std::cout << "path_mm" << "\t" << p1 << "\t" << offset(p2) << "\t" << (is_rev(p2)?"-":"+") << std::endl; });
    }

    // 4) construct the links of the graph in L by rewriting the forward and reverse of Q in terms of pairs of basis in S
    std::string link_fwd_mm_idx = args::get(base) + ".sqlf";
    std::string link_rev_mm_idx = args::get(base) + ".sqlr";
    std::remove(link_fwd_mm_idx.c_str());
    std::remove(link_rev_mm_idx.c_str());
    dmultimap<pos_t, pos_t> link_fwd_mm(link_fwd_mm_idx);
    dmultimap<pos_t, pos_t> link_rev_mm(link_rev_mm_idx);
    derive_links(seqidx, graph_length, path_mm, link_fwd_mm, link_rev_mm);
    if (args::get(debug)) {
        link_fwd_mm.for_each_pair([&](const pos_t& p1, const pos_t& p2) {
                std::cout << "link_fwd_mm" << "\t"
                          << pos_to_string(p1) << "\t"
                          << pos_to_string(p2) << std::endl; });
        link_rev_mm.for_each_pair([&](const pos_t& p1, const pos_t& p2) {
                std::cout << "link_rev_mm" << "\t"
                          << pos_to_string(p1) << "\t"
                          << pos_to_string(p2) << std::endl; });
    }

    // 5) generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    sdsl::bit_vector seq_id_bv(graph_length);
    compact_nodes(graph_length, link_fwd_mm, link_rev_mm, seq_id_bv);
    if (args::get(debug)) std::cerr << seq_id_bv << std::endl;
    sdsl::sd_vector<> seq_id_cbv;
    sdsl::sd_vector<>::rank_1_type seq_id_cbv_rank;
    sdsl::sd_vector<>::select_1_type seq_id_cbv_select;
    sdsl::util::assign(seq_id_cbv, sdsl::sd_vector<>(seq_id_bv));
    seq_id_bv = sdsl::bit_vector(); // clear bitvector
    sdsl::util::assign(seq_id_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_id_cbv));
    sdsl::util::assign(seq_id_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_id_cbv));

    // 6) emit the graph in GFA or VGP format
    if (!args::get(gfa_out).empty()) {
        std::ofstream out(args::get(gfa_out).c_str());
        emit_gfa(out, graph_length, seq_v_file, path_mm, link_fwd_mm, link_rev_mm, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx);
    } else if (!args::get(vgp_base).empty()) {
        emit_vgp(args::get(vgp_base), graph_length, seq_v_file, path_mm, link_fwd_mm, link_rev_mm, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx);
    } else {
        emit_gfa(std::cout, graph_length, seq_v_file, path_mm, link_fwd_mm, link_rev_mm, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx);
    }

    if (!args::get(keep_temp_files)) {
        seqidx.remove_index_files();
        std::remove(aln_idx.c_str());
        std::remove(seq_v_file.c_str());
        std::remove(node_mm_idx.c_str());
        std::remove(path_mm_idx.c_str());
        std::remove(link_fwd_mm_idx.c_str());
        std::remove(link_rev_mm_idx.c_str());
    }

    return(0);
}
