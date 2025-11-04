# Phase 2: Core Data Structures Migration - Research

## Overview

Phase 1 successfully migrated all simple utility modules (Steps 1-23). Phase 2 requires migrating the core data structures that all remaining C++ code depends on.

## Current C++ Data Structures

### 1. mmmulti::iitree (Memory-Mapped Interval Tree)
**Location:** `deps/mmmulti/src/mmiitree.hpp`

**Purpose:** Efficient interval overlap queries using implicit binary tree in sorted array

**Key Features:**
- Memory-mapped storage for large datasets
- Interval overlap queries: O(log n + k) where k = number of overlaps
- Write-once, query-many pattern
- Uses atomic queue for concurrent writes
- Requires indexing step after all intervals added

**Current Implementation:**
- Stores intervals as sorted array on disk
- Implicit binary tree structure (no pointers)
- Each node stores: start, end, max (for subtree), data
- Uses mio::mmap for memory mapping
- Uses ips4o for parallel sorting

### 2. mmmulti::map (Disk-Backed Multimap)
**Location:** `deps/mmmulti/src/mmmultimap.hpp`

**Purpose:** Multimap with numeric keys, memory-mapped storage

**Key Features:**
- Sorted key/value pairs on disk
- Succinct bitvector marks first occurrence of each key
- O(1) key lookup using select queries
- Supports iteration over values for a key

**Current Implementation:**
- Uses sdsl::sd_vector for succinct bitvector storage
- Uses ips4o for sorting
- Memory-mapped with mio::mmap
- Atomic queue for concurrent appends

### 3. mmmulti::set (Disk-Backed Multiset)
**Location:** `deps/mmmulti/src/mmmultiset.hpp`

**Purpose:** Sorted multiset with duplicate values, disk-backed

**Key Features:**
- Sorted values on disk
- Supports iteration and counting

**Current Implementation:**
- Similar to mmmulti::map but simpler (no keys)
- Uses ips4o for sorting
- Memory-mapped with mio::mmap

### 4. seqindex_t (Sequence Index)
**Location:** `src/seqindex.hpp`, `src/seqindex.cpp`

**Purpose:** Index FASTA/FASTQ files for random access to sequences

**Key Features:**
- Parses FASTA/FASTQ
- Concatenates all sequences into single file
- Provides random access by sequence name or position
- Tracks sequence boundaries
- Reverse complement support via pos_t encoding

**Current Implementation:**
- Uses sdsl::sd_vector to mark sequence start positions
- Uses sdsl::csa_wt (compressed suffix array) to index sequence names
- Memory-mapped sequence data
- ~313 lines of code

## Rust Crate Research

### Interval Trees

#### ✅ **rust-lapper** (RECOMMENDED for performance)
- **URL:** https://github.com/sstadick/rust-lapper
- **Performance:** 4-10x faster than other methods
- **Algorithm:** BITS algorithm
- **Features:** serde support, optimized for genomics
- **Cons:** In-memory only (would need custom serialization)

#### **rust-bio interval tree**
- **URL:** https://github.com/rust-bio/rust-bio
- **Implementation:** AVL tree based
- **Features:** Part of comprehensive bioinformatics library
- **Cons:** May not support memory-mapped storage directly

#### **store-interval-tree**
- **URL:** https://crates.io/crates/store-interval-tree
- **Features:** Balanced tree, supports open/closed/unbounded intervals
- **Based on:** rudac and bio

### Succinct Data Structures

#### ✅ **vers-vecs** (RECOMMENDED)
- **URL:** https://crates.io/crates/vers-vecs
- **Performance:** Among fastest publicly available rank/select implementations
- **Features:**
  - Succinct bit vectors with rank/select
  - Pure Rust
  - BMI2 and popcnt support for 2-3x speedup on x86_64
  - Well-maintained

#### **simple-sds**
- **URL:** https://github.com/jltsiren/simple-sds
- **Author:** Jouni Sirén (author of GBWT, expert in succinct structures)
- **Features:**
  - Plain BitVector with rank/select
  - rank(), rank_zero(), select(), select_zero()
  - predecessor(), successor()
  - Designed for bioinformatics

#### **succinct** (deprecated)
- **URL:** https://github.com/tov/succinct-rs
- **Status:** Deprecated but still functional
- **Features:** O(lg lg n) select via binary search over ranks

### Memory Mapping

#### ✅ **memmap2** (STANDARD)
- **URL:** https://crates.io/crates/memmap2
- **Downloads:** 155+ million
- **Status:** Active fork of memmap-rs
- **Features:**
  - Cross-platform
  - Mmap and MmapMut for read/write
  - Safe abstractions

### Sorting

#### ✅ **rayon** (STANDARD)
- **URL:** https://crates.io/crates/rayon
- **Purpose:** Data-parallelism library
- **Features:**
  - par_sort() for parallel sorting
  - Work-stealing scheduler
  - Easy to use

#### **pdqsort**
- Pattern-defeating quicksort (used by Rust std)
- Not explicitly parallel but very fast

### Suffix Arrays & Compressed Indexes

#### **rust-bio suffix arrays**
- Part of rust-bio
- Not compressed like sdsl::csa_wt
- Good for smaller datasets

#### ⚠️ **No direct equivalent to sdsl::csa_wt**
- sdsl (Succinct Data Structure Library) has no complete Rust port
- May need to:
  - Keep C++ version via FFI
  - Use simpler hash-based name lookup
  - Implement minimal CSA ourselves
  - Use FM-index from rust-bio

## Migration Strategy Options

### Option A: Hybrid Approach (RECOMMENDED)
Keep complex C++ structures (mmmulti, sdsl) as-is, migrate only the algorithms

**Pros:**
- Lower risk
- Proven performance
- Faster migration

**Cons:**
- Not pure Rust
- Still have C++ dependencies

### Option B: Pure Rust with Alternatives
Replace all structures with Rust equivalents

**Interval trees:** Custom memory-mapped wrapper around rust-lapper
**Succinct vectors:** vers-vecs or simple-sds
**Sequence index:** Simplified version with hash-based name lookup
**Memory mapping:** memmap2

**Pros:**
- Pure Rust (long-term maintainability)
- Memory safety guarantees
- Can optimize for our specific use case

**Cons:**
- More work upfront
- Need to validate performance
- May need to implement missing features

### Option C: Gradual Replacement
Start with hybrid, incrementally replace structures

1. Use C++ structures via FFI initially
2. Profile to find bottlenecks
3. Replace one structure at a time with Rust
4. Validate performance at each step

**Pros:**
- Maintains working system at every step (Ship of Theseus!)
- Can optimize where it matters
- Reduces risk

**Cons:**
- Longest timeline
- Complex FFI boundaries

## Recommended Approach for seqindex_t

Since seqindex_t is the foundation, we should start there:

### Plan for seqindex Migration

**Phase 2.1: Sequence Storage**
- Use memmap2 for sequence data
- Keep FASTA/FASTQ parsing in Rust
- Store concatenated sequences in memory-mapped file

**Phase 2.2: Sequence Boundaries**
- Use vers-vecs for succinct bitvector marking sequence starts
- Implement rank/select for boundary queries
- Test performance vs sdsl::sd_vector

**Phase 2.3: Name Indexing**
- **Option A:** Simple HashMap<String, usize> for name -> id
  - Fast, simple, but uses more memory
- **Option B:** Keep sdsl::csa_wt via FFI temporarily
  - Proven, but keeps C++ dependency
- **Option C:** Implement minimal suffix array
  - Educational, but time-consuming

**Recommended:** Start with Option A (HashMap), profile memory usage, decide if compression needed

## Estimated Effort

### seqindex_t Migration
- **Time:** 1-2 sessions
- **Complexity:** Medium
- **Risk:** Low (well-defined interface)
- **Lines of code:** ~300-400 Rust (similar to C++)

### mmmulti structures
- **Time:** 3-5 sessions (if building custom)
- **Complexity:** High
- **Risk:** Medium (performance critical)
- **Lines of code:** ~800-1200 Rust

### Total Phase 2
- **Time:** 5-10 sessions
- **Complexity:** High
- **Dependencies:** Need to research and integrate multiple crates

## Next Steps

1. ✅ Complete this research document
2. Start seqindex_t migration with simplified approach:
   - memmap2 for sequence storage
   - vers-vecs for boundary marking
   - HashMap for name lookup
3. Build tests comparing against C++ version
4. Profile memory and performance
5. Decide on mmmulti strategy based on seqindex results

## Open Questions

1. **Memory overhead:** How much memory does HashMap use vs CSA for sequence names?
2. **Performance:** Is vers-vecs rank/select as fast as sdsl::sd_vector?
3. **Disk I/O:** Can we match mmmulti's memory-mapped performance?
4. **Concurrency:** Do we need the atomic queue pattern or can we simplify?

## Crate Dependencies to Add

```toml
[dependencies]
memmap2 = "0.9"           # Memory mapping
vers-vecs = "1.0"         # Succinct bitvectors with rank/select
rayon = "1.10"            # Parallel sorting and iteration
# rust-lapper = "1.1"     # (Maybe later for interval trees)
# bio = "2.0"             # (Alternative: rust-bio)
```

## References

- rust-lapper: https://github.com/sstadick/rust-lapper
- vers-vecs: https://crates.io/crates/vers-vecs
- simple-sds: https://github.com/jltsiren/simple-sds
- memmap2: https://crates.io/crates/memmap2
- Rust-Bio: https://rust-bio.github.io/
- sdsl (C++): https://github.com/simongog/sdsl-lite
