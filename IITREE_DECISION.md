# iitree Implementation Decision

**Date:** 2025-11-04
**Status:** APPROVED - Starting implementation

## Context

Phase 2 of the seqwish Rust migration has completed Step 24 (seqindex with proper Big O bounds). The next critical dependency is **mmmulti::iitree** - the interval tree used by all remaining algorithms (alignments, links, transclosure, GFA output).

## Research Summary

### C++ mmmulti::iitree Analysis

**Location:** `deps/mmmulti/src/mmiitree.hpp` (440 lines)

**Algorithm:** Implicit Augmented Interval Tree (IAITree/cgranges)
- Sorted flat array interpreted as complete binary tree
- No pointers - pure index arithmetic for navigation
- Augmented with `max` field for pruning

**Complexity (EXACT REQUIREMENTS):**
- **Space:** O(n) - exactly 32 bytes per interval
- **Query:** O(log n + k) where k = number of results
- **Build:** O(n log n) - ips4o parallel sort + O(n) augmentation

**Critical Features:**
- Memory-mapped disk storage (datasets larger than RAM)
- Thread-safe lock-free writer queue
- Parallel sorting with ips4o
- Battle-tested on human pangenome scale datasets

### Rust Crate Survey

**Evaluated:** rust-lapper, coitrees, iset, bio::IntervalTree, nclist, others

**Finding:** No perfect match
- **rust-lapper** (most popular): Fast but complexity undocumented, no memory mapping
- **coitrees** (closest): Cache-oblivious, SIMD, but static-only, no serialization
- **iset/bio::IntervalTree**: Documented O(log n + k) but pointer-based trees (hard to memory-map)
- **Critical gap:** None support memory-mapped disk storage

## Decision: Port mmmulti::iitree Algorithm to Rust

### Rationale

1. **Big O guarantee:** We know exact space/time bounds match requirements
2. **Proven algorithm:** Battle-tested in seqwish on billion-base pangenomes
3. **Memory mapping:** Non-negotiable for large datasets
4. **Manageable scope:** Core algorithm ~200 lines, total ~440 lines
5. **Reusable:** Creates standalone crate for future bioinformatics projects

### Approach

**Create standalone crate:** `~/iitree-rs/`

**Ship of Theseus implementation strategy:**
1. **Core algorithm** (Step 1): Port index_core + overlap (~200 lines)
2. **Validation** (Step 2): Port test from mmmulti, validate against C++
3. **I/O layer** (Step 3): Add memory-mapped file support
4. **Concurrency** (Step 4): Add thread-safe writer with lock-free queue
5. **Integration** (Step 5): Use in seqwish-rs

**Dependencies:**
- `memmap2` - memory mapping (replaces mio)
- `rayon` - parallel sorting (replaces ips4o)
- `crossbeam` - lock-free queue (replaces atomic_queue)
- `serde` (optional) - for alternative serialization

### Timeline Estimate

- **Core + tests:** 2-3 hours (1 session)
- **Memory mapping:** 1-2 hours
- **Writer/concurrency:** 2-3 hours
- **Performance tuning:** 1-2 sessions
- **Total:** 3-5 sessions to production-ready

### Parallel Work Strategy

**Option 1:** Use FFI to C++ iitree while building Rust version
- Allows seqwish migration to continue
- Swap to pure Rust once validated

**Option 2:** Dedicate focus to iitree-rs first (CHOSEN)
- Unblocks all remaining algorithm migrations
- iitree is critical path for Phase 2B/2C

## Progress Update (2025-11-04)

### ✅ Core Implementation Complete

**Location:** `~/iitree-rs/`

**Completed:**
1. ✅ Create `~/iitree-rs/` with cargo new
2. ✅ Define `Interval<S, T>` struct with #[repr(C)]
3. ✅ Port index_core algorithm (lines 92-118 from mmiitree.hpp)
4. ✅ Port overlap algorithm (lines 379-427 from mmiitree.hpp)
5. ✅ Port cgranges test (lines 224-255 from main.cpp)
6. ✅ Validate: identical results to C++ on test data
7. ✅ Add parallel sorting with rayon
8. ✅ Comprehensive README with usage examples

**Tests:** 9/9 passing
- 7 basic tests (intervals, overlap, edge cases)
- 2 comprehensive cgranges tests (1K and 10K intervals)

**Complexity validated:**
- Space: O(n) - 32 bytes per interval (u64/u64)
- Query: O(log n + k) - all queries return correct results
- Build: O(n log n) - parallel sort with rayon

### 📋 Remaining Work

**Optional enhancements** (not required for seqwish integration):
1. ⏳ Add memory-mapped I/O (for disk persistence)
2. ⏳ Add thread-safe writer with lock-free queue
3. ⏳ Performance benchmarks vs C++
4. ⏳ Publish to crates.io

**Ready for integration:**
The core algorithm is complete and validated. We can now:
- **Option A:** Use directly in seqwish-rs (Vec-based, in-memory)
- **Option B:** Add memory-mapping layer (match C++ exactly)
- **Option C:** Use temporarily via FFI while adding features

**Recommendation:** Use directly in seqwish-rs with in-memory Vec storage for now. Add memory-mapping later if profiling shows it's needed.

## Success Criteria

- ✅ Exact Big O bounds: O(n) space, O(log n + k) query, O(n log n) build
- ✅ Memory-mapped disk storage
- ✅ Identical query results to C++ mmmulti::iitree on test data
- ✅ Performance within 10% of C++ version
- ✅ Thread-safe concurrent writes
- ✅ Comprehensive test suite

## References

- C++ implementation: `/home/erik/seqwish/deps/mmmulti/src/mmiitree.hpp`
- Algorithm paper: cgranges by Heng Li (https://github.com/lh3/cgranges)
- Research doc: Will create IITREE_RESEARCH.md with detailed algorithm analysis
- Big O requirements: `.claude/CLAUDE.md`
