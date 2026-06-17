#include <iostream>
#include <cstdint>
#include <cstdio>
#include <string>
#include <chrono>
#include <iomanip>
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
#include "exists.hpp"
#include "time.hpp"
#include "utils.hpp"
#include "version.hpp"
#include "tempfile.hpp"

using namespace seqwish;

int main(int argc, char** argv) {
    args::ArgumentParser parser("seqwish: a variation graph inducer\n" + seqwish::Version::get_version() + ": " + seqwish::Version::get_codename());
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> paf_alns(parser, "FILE", "Induce the graph from these PAF formatted alignments. Optionally, a list of filenames and minimum match lengths: [file_1][:min_match_length_1],... This allows the differential filtering of short matches from some but not all inputs, in effect allowing `-k` to be specified differently for each input.", {'p', "paf-alns"});
    args::ValueFlag<std::string> seqs(parser, "FILE", "The sequences used to generate the alignments (FASTA, FASTQ, .seq)", {'s', "seqs"});
    args::ValueFlag<std::string> tmp_base(parser, "PATH", "directory for temporary files [default: `pwd`]", {'b', "temp-dir"});
    args::ValueFlag<std::string> gfa_out(parser, "FILE", "Write the graph in GFA to FILE", {'g', "gfa"});
    args::ValueFlag<std::string> sml_in(parser, "FILE", "Use the sequence match list in FILE to subset the input alignments", {'m', "match-list"});
    //args::ValueFlag<std::string> vgp_base(parser, "BASE", "Write the graph in VGP format with basename FILE", {'o', "vgp-out"});
    args::ValueFlag<int> thread_count(parser, "N", "Use this many threads during parallel steps", {'t', "threads"});
    args::ValueFlag<uint64_t> repeat_max(parser, "N", "Limit transitive closure to include no more than N copies of a given input base", {'r', "repeat-max"});
    args::ValueFlag<uint64_t> min_repeat_dist(parser, "N", "Prevent transitive closure for bases at least this far apart in input sequences", {'l', "min-repeat-distance"});
    args::ValueFlag<uint64_t> min_match_len(parser, "N", "Filter exact matches below this length. This can smooth the graph locally and prevent the formation of complex local graph topologies from forming due to differential alignments.", {'k', "min-match-len"});
    args::ValueFlag<float> match_sparsification(parser, "N", "Sparsify input matches, keeping the fraction that minimize a hash function.", {'f', "sparse-factor"});
    args::ValueFlag<std::string> transclose_batch(parser, "N", "Number of bp to use for transitive closure batch (1k = 1K = 1000, 1m = 1M = 10^6, 1g = 1G = 10^9) [default 1M]", {'B', "transclose-batch"});
    //args::ValueFlag<uint64_t> num_domains(parser, "N", "number of domains for iitii interpolation", {'D', "domains"});
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during graph induction", {'T', "keep-temp"});
    args::Flag show_progress(parser, "show-progress", "log algorithm progress", {'P', "show-progress"});
    args::Flag verbose_debug(parser, "verbose-debug", "enable verbose debugging", {'V', "verbose-debug"});
	args::Flag version(parser, "version", "show the current version including github commit hash", {'v', "version"});
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

	if (version) {
		std::cerr << seqwish::Version::get_version() << std::endl;
		exit(0);
	}

    std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
    
    size_t num_threads = args::get(thread_count) ? args::get(thread_count) : 1;
    /*
    if (num_threads) {
        omp_set_num_threads(args::get(thread_count));
    } else {
        omp_set_num_threads(1);
    }
    */

    if (!args::get(seqs).empty() && !file_exists(args::get(seqs))) {
        std::cerr << "[seqwish] ERROR: input sequence file " << args::get(seqs) << " does not exist" << std::endl;
        return 2;
    }

    // parse paf args
    std::vector<std::pair<std::string, uint64_t>> pafs_and_min_lengths;
    if (!args::get(paf_alns).empty()) {
        pafs_and_min_lengths = parse_paf_spec(args::get(paf_alns));
        for (auto& p : pafs_and_min_lengths) {
            if (!file_exists(p.first)) {
                std::cerr << "[seqwish] ERROR: input alignment file " << args::get(paf_alns) << " does not exist" << std::endl;
                return 4;
            }else {
                 // Check if the first non-empty line has the CIGAR

                igzstream paf_in(p.first.c_str());

                while (!paf_in.eof()) {
                    std::string line;
                    std::getline(paf_in, line);

                    if (!line.empty()) {
                        paf_row_t paf(line);

                        if (paf.cigar.empty()){
                            std::cerr << "[seqwish] WARNING: input alignment file " << p.first << " does not have CIGAR strings. "
                                      << "The resulting graph will only represent the input sequences." << std::endl;
                        }
                        break;
                    }
                }
            }
        }
    }

    if (tmp_base) {
        temp_file::set_dir(args::get(tmp_base));
    } else {
        char cwd[512];
        getcwd(cwd, sizeof(cwd));
        temp_file::set_dir(std::string(cwd));
    }

    temp_file::set_keep_temp(args::get(keep_temp_files));

    // 1) index the queries (Q) to provide sequence name to position and position to sequence name mapping, generating a CSA and a sequence file
    if (args::get(show_progress)) std::cerr << "[seqwish::seqidx] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " indexing sequences" << std::endl;
    auto seqidx_ptr = std::make_unique<seqindex_t>();
    auto& seqidx = *seqidx_ptr;
    seqidx.build_index(args::get(seqs));
    seqidx.save();
    if (args::get(show_progress)) std::cerr << "[seqwish::seqidx] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " index built" << std::endl;

    // 2) parse the alignments into position pairs and index (A)
    if (args::get(show_progress)) std::cerr << "[seqwish::alignments] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " processing alignments" << std::endl;
    const std::string aln_idx = temp_file::create("seqwish-", ".sqa");
    auto aln_iitree_ptr = std::make_unique<mmmulti::iitree<uint64_t, pos_t>>(aln_idx);
    auto& aln_iitree = *aln_iitree_ptr;
    aln_iitree.open_writer();
    float sparse_match = match_sparsification ? args::get(match_sparsification) : 0;
    if (!pafs_and_min_lengths.empty()) {
        for (auto& p : pafs_and_min_lengths) {
            auto& file = p.first;
            uint64_t min_length = p.second;
            if (!min_length && args::get(min_match_len)) {
                min_length = args::get(min_match_len);
            }
            unpack_paf_alignments(file, aln_iitree, seqidx, min_length, sparse_match, num_threads);
        }
    }
    if (args::get(show_progress)) std::cerr << "[seqwish::alignments] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " indexing" << std::endl;
    aln_iitree.index(num_threads);
    if (args::get(show_progress)) std::cerr << "[seqwish::alignments] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " index built" << std::endl;
    //if (args::get(debug)) dump_paf_alignments(args::get(paf_alns));
    //uint64_t n_domains = std::max((uint64_t)1, (uint64_t)args::get(num_domains));
    //range_pos_iitii aln_iitree = aln_iitree_builder.build(n_domains);

    if (args::get(verbose_debug)) {
        for (auto& interval : aln_iitree) {
            std::cerr << "aln_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
    }

    // 3) find the transitive closures via the alignments and construct the graph sequence S, and the N and P interval sets
    const std::string seq_v_file = temp_file::create("seqwish-", ".sqs");
    const std::string node_iitree_idx = temp_file::create("seqwish-", ".sqn");
    const std::string path_iitree_idx = temp_file::create("seqwish-", ".sqp");
    auto node_iitree_ptr = std::make_unique<mmmulti::iitree<uint64_t, pos_t>>(node_iitree_idx); // maps graph seq to input seq
    auto& node_iitree = *node_iitree_ptr;
    auto path_iitree_ptr = std::make_unique<mmmulti::iitree<uint64_t, pos_t>>(path_iitree_idx); // maps input seq to graph seq
    auto& path_iitree = *path_iitree_ptr;
    if (args::get(show_progress)) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " computing transitive closures" << std::endl;
    size_t graph_length = compute_transitive_closures(seqidx, aln_iitree, seq_v_file, node_iitree, path_iitree,
                                                      args::get(repeat_max),
                                                      args::get(min_repeat_dist),
                                                      transclose_batch ? (uint64_t)seqwish::handy_parameter(args::get(transclose_batch), 1000000) : 1000000,
                                                      args::get(show_progress),
                                                      num_threads,
                                                      start_time);
    if (args::get(show_progress)) std::cerr << "[seqwish::transclosure] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " done with transitive closures" << std::endl;

    if (args::get(verbose_debug)) {
        for (auto& interval : node_iitree) {
            std::cerr << "node_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
        for (auto& interval : path_iitree) {
            std::cerr << "path_iitree " << interval.st << "-" << interval.en << " " << pos_to_string(interval.data) << std::endl;
        }
    }

    // 4) generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    if (args::get(show_progress)) std::cerr << "[seqwish::compact] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " compacting nodes" << std::endl;
    sdsl::bit_vector seq_id_bv(graph_length+1);
    compact_nodes(seqidx, graph_length, node_iitree, path_iitree, seq_id_bv, num_threads);
    if (args::get(show_progress)) std::cerr << "[seqwish::compact] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " done compacting" << std::endl;
    if (args::get(verbose_debug)) std::cerr << seq_id_bv << std::endl;
    sdsl::sd_vector<> seq_id_cbv;
    sdsl::sd_vector<>::rank_1_type seq_id_cbv_rank;
    sdsl::sd_vector<>::select_1_type seq_id_cbv_select;
    sdsl::util::assign(seq_id_cbv, sdsl::sd_vector<>(seq_id_bv));
    seq_id_bv = sdsl::bit_vector(); // clear bitvector
    sdsl::util::assign(seq_id_cbv_rank, sdsl::sd_vector<>::rank_1_type(&seq_id_cbv));
    sdsl::util::assign(seq_id_cbv_select, sdsl::sd_vector<>::select_1_type(&seq_id_cbv));
    if (args::get(show_progress)) std::cerr << "[seqwish::compact] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " built node index" << std::endl;

    // 5) determine links between nodes
    if (args::get(show_progress)) std::cerr << "[seqwish::links] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " finding graph links" << std::endl;
    const std::string link_mm_idx =  temp_file::create("seqwish-", ".sql");
    auto link_mmset_ptr = std::make_unique<mmmulti::set<std::pair<pos_t, pos_t>>>(link_mm_idx);
    auto& link_mmset = *link_mmset_ptr;
    derive_links(seqidx, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, link_mmset, num_threads);
    if (args::get(show_progress)) std::cerr << "[seqwish::links] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " links derived" << std::endl;

    // 6) emit the graph in GFA or VGP format
    if (args::get(show_progress)) std::cerr << "[seqwish::gfa] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " writing graph" << std::endl;
    if (!args::get(gfa_out).empty()) {
        std::ofstream out(args::get(gfa_out).c_str());
        emit_gfa(out, graph_length, seq_v_file, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx, link_mmset, num_threads);
    /*} else if (!args::get(vgp_base).empty()) {
        assert(false);
        //emit_vgp(args::get(vgp_base), graph_length, seq_v_file, path_mm, link_fwd_mm, link_rev_mm, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx);
    */} else {
        emit_gfa(std::cout, graph_length, seq_v_file, node_iitree, path_iitree, seq_id_cbv, seq_id_cbv_rank, seq_id_cbv_select, seqidx, link_mmset, num_threads);
    }
    if (args::get(show_progress)) std::cerr << "[seqwish::gfa] " << std::fixed << std::showpoint << std::setprecision(3) << seconds_since(start_time) << " done" << std::endl;

    return(0);
}
