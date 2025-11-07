# seqwish 🦀

[![Build Status](https://github.com/pangenome/seqwish/workflows/build%20and%20test/badge.svg)](https://github.com/pangenome/seqwish/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A variation graph inducer** - Build pangenome graphs from pairwise sequence alignments.

Seqwish implements a lossless conversion from pairwise alignments between sequences to a variation graph encoding the sequences and their alignments. As input we typically take all-versus-all alignments, but the exact structure of the alignment set may be defined in an application specific way.

## ✨ Features

- **Memory-safe** - Rewritten in Rust with compile-time safety guarantees
- **Parallel** - Multi-threaded processing throughout the pipeline
- **Scalable** - Disk-backed data structures for processing large genomes
- **Fast** - Performance comparable to the highly-optimized C++ version
- **Verified** - Produces byte-for-byte identical output to reference implementation

## 🚀 Installation

```bash
# From source (recommended)
git clone https://github.com/pangenome/seqwish
cd seqwish
cargo build --release
# Binary will be in target/release/seqwish

# Future: cargo install seqwish (once published to crates.io)
```

## 📖 Usage

```bash
# Basic usage
seqwish -s sequences.fa -p alignments.paf -g output.gfa

# With options for large datasets
seqwish \
  -s sequences.fa \      # Input sequences (FASTA/FASTQ)
  -p alignments.paf \    # Pairwise alignments (PAF format)
  -g output.gfa \        # Output variation graph (GFA format)
  -t 16 \                # Use 16 threads
  -k 19 \                # Filter matches < 19bp
  -P                     # Show progress
```

### Options

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

## 🔬 Algorithm Overview

The algorithm proceeds in stages:

1. **Sequence Indexing** - Build FM-index of input sequences
2. **Alignment Indexing** - Parse PAF alignments into interval trees
3. **Transitive Closure** - Compute equivalence classes of aligned positions
4. **Graph Emission** - Write graph sequence from closures
5. **Node Compaction** - Merge non-bifurcating regions
6. **Link Derivation** - Extract edges between nodes
7. **GFA Output** - Emit final variation graph

## 📊 Example Workflow

```bash
# 1. Generate all-to-all alignments
minimap2 -cx asm20 -X sequences.fa sequences.fa > alignments.paf

# 2. Build the variation graph
seqwish -s sequences.fa -p alignments.paf -g graph.gfa -t 16 -P

# 3. Visualize (requires vg and graphviz)
vg view -dp graph.gfa | dot -Tpng > graph.png
```

## 🏗️ Use Cases

- **Pangenome construction** - Build graphs from multiple related genomes
- **Structural variation** - Capture large-scale genomic rearrangements
- **Population genomics** - Represent variation across many samples
- **Reference graphs** - Create enhanced reference structures

## 🔗 Related Tools

- [minimap2](https://github.com/lh3/minimap2) - Generate PAF alignments
- [vg](https://github.com/vgteam/vg) - Variation graph toolkit
- [odgi](https://github.com/pangenome/odgi) - Graph optimization and manipulation
- [gfaffix](https://github.com/marschall-lab/GFAffix) - GFA graph manipulation

## 📚 Citation

If you use seqwish, please cite:

```
Garrison E, Guarracino A. (2023)
Unbiased pangenome graphs
Bioinformatics, Volume 39, Issue 1, btac743
https://doi.org/10.1093/bioinformatics/btac743
```

## 🏛️ Implementation

This is a complete Rust reimplementation of seqwish using the "Ship of Theseus" pattern - incrementally replacing components while maintaining correctness. The Rust version:

- ✅ Produces **byte-for-byte identical** output to C++
- ✅ Passes all original test suites
- ✅ Provides memory safety guarantees
- ✅ Offers comparable performance

The original C++ implementation is preserved in `cpp/` for reference. See [README_CPP.md](README_CPP.md) for C++ documentation.

## 🤝 Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

MIT - see [LICENSE](LICENSE) file

## 💬 Support

- 🐛 [Issue Tracker](https://github.com/pangenome/seqwish/issues)
- 💬 [Discussions](https://github.com/pangenome/seqwish/discussions)

---

**Made with 🦀 Rust**
