# seqwish

[![build and test](https://github.com/pangenome/seqwish/actions/workflows/build_and_test_on_push.yml/badge.svg)](https://github.com/pangenome/seqwish/actions/workflows/build_and_test_on_push.yml)
[![crates.io](https://img.shields.io/crates/v/seqwish.svg)](https://crates.io/crates/seqwish)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/seqwish/README.html)

*These <b>seq</b>uences <b>wish</b> they were squished into a graph.*

## a variation graph inducer

`seqwish` implements a lossless conversion from pairwise alignments between sequences to a variation graph encoding the sequences and their alignments.
As input we typically take all-versus-all alignments, but the exact structure of the alignment set may be defined in an application specific way.
This algorithm uses a series of disk-backed sorts and passes over the alignment and sequence inputs to allow the graph to be constructed from very large inputs that are commonly encountered when working with large numbers of noisy input sequences.
Memory usage during construction and traversal is limited by the use of sorted disk-backed arrays and succinct rank/select dictionaries to record a queryable version of the graph.

## Rust implementation

This version of seqwish has been completely rewritten in Rust. The new implementation:

- Produces **byte-for-byte identical output** to the original C++ version
- Provides a **library API** for programmatic use (`use seqwish;`)
- Supports both **disk-backed** and **in-memory** operation modes
- Is available on [crates.io](https://crates.io/crates/seqwish)

The original C++ implementation is preserved in `cpp/` for reference.

## Citation

Erik Garrison, Andrea Guarracino, **Unbiased pangenome graphs**, Bioinformatics, Volume 39, Issue 1, January 2023, btac743, https://doi.org/10.1093/bioinformatics/btac743

## squish graph induction algorithm

As input we have *Q*, which is a concatenation of the sequences from which we will build the graph.
We build an compressed suffix array (CSA) mapping sequence names to offsets in *Q*, and also the inverse using a rank/select dictionary on a bitvector marking the starts of sequences in *Q*.
This allows us to map between positions in the sequences of *Q*, which is the format in which alignment algorithms typically express alignments, and positions in *Q* itself, which is the coordinate space we will use as a basis for the generation of our graph.
To relate the sequences in *Q* to each other we apply a function *map* to generate alignments *A*.
Although these alignments tend to be represented using oriented interval pairs in *Q*, for simplicity and robustness to graph complexity, we describe *A* as a vector of pairs of bidirectional positions (sequence offsets and strands) *b* in *Q* , such that *A* = [(*b<sub>q</sub>*, *b<sub>r</sub>*), ... ].
We sort *A* by the first member (*b<sub>q</sub>*) of each pair, ensuring that the entries in *A* are ordered according to their order in *Q*.

To query the induced graph we build a rank/select dictionary allowing efficient traversal of *A*, based on a bit vector *A<sub>bv</sub>* of the same length as *A* such that we record a 1 at those positions which correspond to the first instance of a given *b<sub>q</sub>* and record a 0 in *A<sub>bv</sub>* otherwise.
We record which *b<sub>q</sub>* we have processed in the bitvector *Q<sub>seen</sub>* which is of the same length as *Q*.
This allows us to avoid a quadratic penalty in the order of the size of the transitive closures in *Q* given by the *map* function.

Now we inductively derive the graph implied by the alignments.
For each base *b<sub>q</sub>* in *Q*, we find its transitive closure *c<sub>q</sub>* := [*b<sub>q</sub>*, *b<sub>r<sub>1</sub></sub>*, ... ] via the *map* operation by traversing the aligned base pairs recorded in *A*.
We write the character of the base *b<sub>q</sub>* to a vector *S*, then for each *b<sub>c</sub>* in *c<sub>q</sub>* we record a pair [*s<sub>i</sub>, b<sub>c</sub>*] into *N* and its reverse, [*b<sub>c</sub>*, *s<sub>i</sub>*] into *P*.
We mark *Q<sub>seen</sub>* for each base in each emitted cluster, and we do not consider marked bases in subsequent transitive closures.
By sorting *N* and *P* by their first entries, we can build rank/select dictionaries on them akin to that we built on *A* that allow random access by graph base (as given in *S*) or input base (as given in *Q*).

To fully induce the variation graph we need to establish the links between bases in *S* that would be required for us to find any sequence in the input as a walk through the graph.
We do so by rewriting *Q* (in both the forward and reverse orientation) in terms of pairs of bases in *S*, then sorting the resulting pairs by their first element, which yields *L* = [(*b<sub>a</sub>*, *b<sub>b</sub>*), ... ].
These pairs record the links and their frequencies, which we can emit or filter (such as by frequency) as needed given particular applications.
In typical use we take the graph to be given by the unique elements of *L*.

Our data model encodes the graph using single-base nodes, but often downstream use requires identifying nodes and thus we benefit from compressing the unitigs of the graph into single nodes, which reduces memory used by identifiers in analysis.
We can compress the node space of the graph by traversing *S*, and for each base querying the inbound links.
Maintaining a bitvector *S<sub>id</sub>* of length equal to *S* we mark each base at which we see any link other than one from or to the previous base on the forward or reverse strand, or at bases where we have no incoming links.
By building a rank/select dictionary on *S<sub>id</sub>* we can assign a smaller set of node ids to the sequence space of the graph.

Given the id space encoded by *S<sub>id</sub>* we can materialize the graph in a variety of interchange formats, or provide id-based interfaces to the indexed squish graph.
To generate graphs in [.vg](https://github.com/vgteam/vg/blob/master/src/vg.proto) or [GFAv1 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), we want to decompose the graph into its nodes (*S*), edges (*L*) and paths (*P*).
The nodes are given by *S* and *S<sub>id</sub>*, and similarly we project *L* and *P* through *S<sub>id</sub>* to obtain a compressed variation graph.

## using our graph

At present we rest on existing tools, such as [vg](https://github.com/vgteam/vg), [GCSA2](https://github.com/jltsiren/gcsa2), [GBWT](https://github.com/jltsiren/gbwt), and other GFAv1-compatible tools to support downstream uses for the graph.
We observe that it will be trivial to use this data model as the basis for various graph traversal and modification algorithms.
Implementation of these features will follow.
We plan to link into the [handle graph model](https://github.com/vgteam/vg/blob/master/src/handle.hpp) developed in the vg project, which should allow the application of any algorithms defined using that interface on the backing graph developed here.
Users familiar with concepts in assembly graphs will notice many similarities between the graph induction algorithm and typical steps in assembly algorithms, and we intend to explore these as well using this platform.

## building

### from source (Rust)

You'll need Rust installed. Then:

```bash
git clone https://github.com/pangenome/seqwish.git
cd seqwish
cargo build --release
# Binary will be in target/release/seqwish
```

### from crates.io

```bash
cargo install seqwish
```

### Bioconda

`seqwish` recipes for Bioconda are available at https://bioconda.github.io/recipes/seqwish/README.html.
To install the latest version using `Conda` execute:

```bash
conda install -c bioconda seqwish
```

## usage

`seqwish` supports PAF format output of several sequence aligners, like [wfmash](https://github.com/ekg/wfmash) and [minimap2](https://github.com/lh3/minimap2).
It requires the CIGAR string of the alignment to be provided in the `cg:z:` optional field. It uses large temporary files during the construction.
By default, these are prefixed with the output GFA file name, but this can be changed with the `-b[base], --base=[base]` command line argument. The input sequences can be in FASTA or FASTQ format, either in plain text or gzipped.
It writes [GFA1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#the-gfa-format-specification) on its standard output.

### basic usage

```bash
seqwish -s sequences.fa -p alignments.paf -g output.gfa
```

### with options

```bash
seqwish \
  -s sequences.fa \      # Input sequences (FASTA/FASTQ)
  -p alignments.paf \    # Pairwise alignments (PAF format)
  -g output.gfa \        # Output variation graph (GFA format)
  -t 16 \                # Use 16 threads
  -k 19 \                # Filter matches < 19bp
  -P                     # Show progress
```

### options

```
-s, --seqs <FILE>              Input sequences (FASTA/FASTQ, optionally gzipped)
-p, --paf-alns <FILE>          Input alignments (PAF format, optionally gzipped)
-g, --gfa <FILE>               Output graph (GFA v1.0 format)
-t, --threads <N>              Number of threads [default: 1]
-k, --min-match-len <N>        Minimum match length [default: 0]
-r, --repeat-max <N>           Maximum repeat copies in transitive closure [default: 0]
-l, --min-repeat-distance <N>  Minimum distance for repeat handling [default: 0]
-B, --transclose-batch <N>     Transitive closure batch size [default: 1000000]
-b, --temp-dir <PATH>          Temporary file directory
-T, --keep-temp                Keep temporary files
-P, --show-progress            Show progress messages
```

### wfmash example

```bash
wfmash x.fa x.fa -X > x.paf
seqwish -s x.fa -p x.paf -g x.gfa
```

### minimap2 example

[minimap2](https://github.com/lh3/minimap2) does not emit the CIGAR string in PAF output by default. To do this, specify the `-c` flag:

```bash
minimap2 x.fa x.fa -c -X > x.paf
seqwish -s x.fa -p x.paf -g x.gfa
```

## library usage

The Rust implementation provides a library API:

```rust
use seqwish::{CompactGraph, SeqIndex, paf};

// Build a graph programmatically
let seqidx = SeqIndex::build_from_file("sequences.fa", "workdir")?;
let alignments = paf::parse_paf_file("alignments.paf", &seqidx)?;
// ... process alignments and build graph
```

See the [API documentation](https://docs.rs/seqwish) for details.

## license

MIT
