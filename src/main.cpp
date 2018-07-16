#include <iostream>
#include <cstdint>
#include <string>
#include "args.hxx"
#include "dmultimap.hpp"
#include "seqindex.hpp"
#include "paf.hpp"
#include "alignments.hpp"

using namespace std;
using namespace seqwish;


int main(int argc, char** argv) {
    args::ArgumentParser parser("seqwish: a variation graph inducer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<string> alns(parser, "alns", "induce the graph from these alignments", {'a', "alns"});
    args::ValueFlag<string> seqs(parser, "seqs", "the sequences used to generate the alignments", {'s', "seqs"});
    args::ValueFlag<string> base(parser, "base", "build graph using this basename", {'b', "base"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        cout << parser;
        return 0;
    } catch (args::ParseError e) {
        cerr << e.what() << endl;
        cerr << parser;
        return 1;
    }
    if (argc==1) {
        cout << parser;
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
    dmultimap<int64_t, int64_t> aln_mm(args::get(base) + ".sqa");
    if (args::get(debug)) dump_alignments(args::get(alns));
    unpack_alignments(args::get(alns), aln_mm, seqidx); // yields A index
    // 3) find the transitive closures via the alignments and construct S, N, and P indexed arrays
    // 4) construct the links of the graph in L by rewriting the forward and reverse of Q in terms of pairs of basis in S
    // 5) optionally generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    // 6) emit the graph in GFA

    return(0);
}
