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

using namespace seqwish;

int main(int argc, char** argv) {
    args::ArgumentParser parser("seqwish: a variation graph inducer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> alns(parser, "alns", "induce the graph from these alignments", {'a', "alns"});
    args::ValueFlag<std::string> seqs(parser, "seqs", "the sequences used to generate the alignments", {'s', "seqs"});
    args::ValueFlag<std::string> base(parser, "base", "build graph using this basename", {'b', "base"});
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

    // 1) index the queries (Q) to provide sequence name to position and position to sequence name mapping, generating a CSA and a sequence file
    seqindex_t seqidx;
    seqidx.build_index(args::get(seqs));
    seqidx.save();
    /*
    if (args::get(debug)) {
        seqidx.to_fasta(cout);
    }
    */
    if (args::get(base).empty()
        || args::get(alns).empty()) {
        // nothing to do
        return 0;
    }

    // 2) parse the alignments into position pairs and index (A)
    std::string aln_idx = args::get(base) + ".sqa";
    std::remove(aln_idx.c_str());
    dmultimap<uint64_t, pos_t> aln_mm(aln_idx);
    if (args::get(debug)) dump_alignments(args::get(alns));
    unpack_alignments(args::get(alns), aln_mm, seqidx); // yields array A

    // 3) find the transitive closures via the alignments and construct S, N, and P indexed arrays
    std::string seq_v_file = args::get(base) + ".sqs";
    std::string node_mm_idx = args::get(base) + ".sqn";
    std::string path_mm_idx = args::get(base) + ".sqp";
    std::remove(node_mm_idx.c_str());
    dmultimap<uint64_t, pos_t> node_mm(node_mm_idx);
    std::remove(path_mm_idx.c_str());
    dmultimap<uint64_t, pos_t> path_mm(path_mm_idx);
    sdsl::bit_vector q_seen_bv(seqidx.seq_length());
    compute_transitive_closures(seqidx, aln_mm, q_seen_bv, seq_v_file, node_mm, path_mm);

    // 4) construct the links of the graph in L by rewriting the forward and reverse of Q in terms of pairs of basis in S
    std::string link_mm_idx = args::get(base) + ".sql";
    std::remove(link_mm_idx.c_str());
    dmultimap<pos_t, pos_t> link_mm(link_mm_idx);
    derive_links(seqidx, path_mm, link_mm);
    
    // 5) generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    // 6) emit the graph in GFA

    return(0);
}
