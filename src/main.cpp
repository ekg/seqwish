#include <iostream>
#include <cstdint>
#include <string>
#include "args.hxx"
#include "dmultimap.hpp"
#include "seqindex.hpp"

using namespace std;
using namespace seqwish;

struct arrayX {
    char bin[100];
};

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
    // index the queries (Q) to provide sequence name to position and position to sequence name mapping, generating a CSA and a sequence file
    // parse the alignments into position pairs and index (A)
    // find the transitive closures via the alignments and construct S, N, and P indexed arrays
    // construct the links of the graph in L by rewriting the forward and reverse of Q in terms of pairs of basis in S
    // optionally generate the node id index (I) by compressing non-bifurcating regions of the graph into nodes
    // emit the graph in GFA
    dmultimap<int64_t, arrayX> d;
    SeqIndex idx;
    idx.build_index(args::get(seqs));
    idx.save();
    if (args::get(debug)) {
        idx.to_fasta(cout);
    }
    //cout << args::get(alns) << endl;
    //cout << args::get(seqs) << endl;
    //cout << args::get(base) << endl;
    return(0);
}
