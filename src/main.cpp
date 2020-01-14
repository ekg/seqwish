#include <iostream>
#include <cstdint>
#include <cstdio>
#include <string>
#include "args.hxx"
#include "mmmultimap.hpp"
#include "mmiitree.hpp"
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
#include "match.hpp"
#include "threads.hpp"
#include "exists.hpp"
//#include "iitii_types.hpp"

using namespace seqwish;

int main(int argc, char** argv) {
    args::ArgumentParser parser("seqwish: a variation graph inducer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> paf_alns(parser, "FILE", "induce the graph from these PAF formatted alignments", {'p', "paf-alns"});
    args::ValueFlag<std::string> seqs(parser, "FILE", "the sequences used to generate the alignments (FASTA, FASTQ, .seq)", {'s', "seqs"});
    args::ValueFlag<std::string> base(parser, "BASE", "build graph using this basename", {'b', "base"});
    args::ValueFlag<std::string> gfa_out(parser, "FILE", "write the graph in GFA to FILE", {'g', "gfa"});
    args::ValueFlag<std::string> sml_in(parser, "FILE", "use the sequence match list in FILE to subset the input alignments", {'m', "match-list"});
    args::ValueFlag<std::string> vgp_base(parser, "BASE", "write the graph in VGP format with basename FILE", {'o', "vgp-out"});
    args::ValueFlag<uint64_t> num_threads(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<uint64_t> repeat_max(parser, "N", "limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    args::ValueFlag<uint64_t> min_match_len(parser, "N", "filter exact matches below this length", {'L', "min-match-len"});
    args::ValueFlag<uint64_t> min_transclose_len(parser, "N", "follow transitive closures only through matches >= this", {'k', "min-transclose-len"});
    args::ValueFlag<uint64_t> transclose_batch(parser, "N", "number of bp to use for transitive closure batch (default 1M)", {'B', "transclose-batch"});
    //args::ValueFlag<uint64_t> num_domains(parser, "N", "number of domains for iitii interpolation", {'D', "domains"});
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during graph induction", {'T', "keep-temp"});
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
    if (!args::get(paf_alns).empty() && !file_exists(args::get(paf_alns))) {
        std::cerr << "[seqwish] ERROR: input alignment file " << args::get(paf_alns) << " does not exist" << std::endl;
        return 4;
    }

    std::string work_base = args::get(base);
    if (work_base.empty()) {
        work_base = args::get(gfa_out);
    }

    // 1) index the queries (Q) to provide sequence name to position and position to sequence name mapping, generating a CSA and a sequence file
    seqindex_t seqidx;
    seqidx.build_index(args::get(seqs), work_base);
    seqidx.save();

    // 2) parse the alignments into position pairs and index (A)
    std::string aln_idx = work_base + ".sqa";
    std::remove(aln_idx.c_str());
    mmmulti::iitree<uint64_t, pos_t> aln_iitree(aln_idx);
    if (!args::get(paf_alns).empty()) {
        unpack_paf_alignments(args::get(paf_alns), aln_iitree, seqidx, args::get(min_match_len));
    }
    //if (args::get(debug)) dump_paf_alignments(args::get(paf_alns));
    //uint64_t n_domains = std::max((uint64_t)1, (uint64_t)args::get(num_domains));
    //range_pos_iitii aln_iitree = aln_iitree_builder.build(n_domains);

    if (args::get(debug)) {
        for (auto& interval : aln_iitree) {
            std::cerr << "aln_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
    }

    // 3) find the transitive closures via the alignments and construct the graph sequence S, and the N and P interval sets
    std::string seq_v_file = work_base + ".sqs";
    std::string node_iitree_idx = work_base + ".sqn";
    std::string path_iitree_idx = work_base + ".sqp";
    std::remove(seq_v_file.c_str());
    std::remove(node_iitree_idx.c_str());
    std::remove(path_iitree_idx.c_str());
    mmmulti::iitree<uint64_t, pos_t> node_iitree(node_iitree_idx); // maps graph seq to input seq
    mmmulti::iitree<uint64_t, pos_t> path_iitree(path_iitree_idx); // maps input seq to graph seq
    size_t graph_length = compute_transitive_closures(seqidx, aln_iitree, seq_v_file, node_iitree, path_iitree,
                                                      args::get(repeat_max), args::get(min_transclose_len),
                                                      !args::get(transclose_batch) ? 1000000 : args::get(transclose_batch));

    if (args::get(debug)) {
        for (auto& interval : node_iitree) {
            std::cerr << "node_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
        for (auto& interval : path_iitree) {
            std::cerr << "path_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
    }

    // 4) generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    sdsl::bit_vector seq_id_bv(graph_length+1);
    compact_nodes(seqidx, graph_length, node_iitree, path_iitree, seq_id_bv);
    if (args::get(debug)) std::cerr << seq_id_bv << std::endl;
    sdsl::sd_vector<> seq_id_cbv;
    sdsl::sd_vector<>::rank_1_type seq_id_cbv_rank;
    sdsl::sd_vector<>::select_1_type seq_id_cbv_select;
    sdsl::util::assign(seq_id_cbv, sdsl::sd_vector<>(seq_id_bv));
    seq_id_bv = sdsl::bit_vector(); // clear bitvector
    sdsl::util::assign(seq_id_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_id_cbv));
    sdsl::util::assign(seq_id_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_id_cbv));

    // 5) determine links between nodes
    std::string link_mm_idx = work_base + ".sql";
    std::remove(link_mm_idx.c_str());
    mmmulti::set<std::pair<pos_t, pos_t>> link_mmset(link_mm_idx);
    derive_links(seqidx, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, link_mmset);
    
    // 6) emit the graph in GFA or VGP format
    if (!args::get(gfa_out).empty()) {
        std::ofstream out(args::get(gfa_out).c_str());
        emit_gfa(out, graph_length, seq_v_file, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx, link_mmset);
    } else if (!args::get(vgp_base).empty()) {
        assert(false);
        //emit_vgp(args::get(vgp_base), graph_length, seq_v_file, path_mm, link_fwd_mm, link_rev_mm, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx);
    } else {
        emit_gfa(std::cout, graph_length, seq_v_file, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx, link_mmset);
    }

    if (!args::get(keep_temp_files)) {
        seqidx.remove_index_files();
        std::remove(aln_idx.c_str());
        std::remove(seq_v_file.c_str());
        std::remove(node_iitree_idx.c_str());
        std::remove(path_iitree_idx.c_str());
        link_mmset.close_reader();
        std::remove(link_mm_idx.c_str());
    }

    return(0);
}
