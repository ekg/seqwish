# Phase 2 Design: Core Data Structures Migration

## Status

**Started:** 2025-11-04
**Current Step:** 24 (seqindex basic implementation complete)

## Overview

Phase 2 focuses on migrating the core data structures that all remaining C++ algorithms depend on. We're using a **hybrid approach**: start with simplified Rust implementations, validate against C++ versions, then incrementally optimize.

## Completed Work

### ✅ Step 24: seqindex - CORRECT Big O Implementation

**Status:** COMPLETE (93 tests passing)

**Implementation:**
- Pure Rust FASTA/FASTQ parser with gzip support
- Memory-mapped sequence storage (memmap2)
- **FM-index/CSA for sequence name lookup** (fm-index 0.3.0)
- **Succinct bitvector for sequence boundaries** (vers-vecs 1.8.1)
- All core query operations implemented
- **Matches C++ sdsl complexity exactly**

**Files:**
- `seqwish-rs/src/seqindex.rs` (574 lines)
- 4 tests covering parsing, access, position queries, and name lookup

**Space Complexity (matching C++ sdsl):**
- Name index: **O(n log σ) bits** using FM-index/CSA (where n = total name chars, σ = alphabet size)
- Sequence boundaries: **O(m log(N/m)) bits** using succinct bitvector (where m = # sequences, N = total length)
- C++ used: sdsl::csa_wt and sdsl::sd_vector
- Rust uses: fm-index::FMIndexWithLocate and vers_vecs::RsVec

**Time Complexity:**
- Name lookup: O(m log n + occ) with FM-index locate()
- Boundary queries: O(1) with succinct select
- Position to seq ID: O(log m) with succinct rank

**Key API Learnings (documented for future use):**
1. fm-index requires text to end with exactly one '\0' character
2. FMIndexWithLocate<C> is generic over character type
3. Search trait must be imported to use iter_matches()
4. MatchWithLocate trait must be imported to use locate()
5. RsVec::rank1(pos) counts 1-bits up to but EXCLUDING pos
6. RsVec::select1(n) returns usize directly, not Option<usize>
7. BitVec uses append_bit() to build, not push()

**Next:** Add FFI bindings so C++ can use it (Step 25)

## Remaining Core Structures

### 📋 mmmulti::iitree (Interval Tree)

**Priority:** HIGH (used by alignments, links, transclosure, GFA output)

**Approach Options:**

#### Option A: Wrap existing C++ via FFI (FAST, keeps proven code)
- Expose mmmulti::iitree as opaque handle
- Add FFI for: add(), overlap(), index(), save/load()
- **Pros:** No porting needed, proven performance
- **Cons:** Still C++ dependency

#### Option B: Custom Rust wrapper around rust-lapper (MEDIUM effort)
- Use rust-lapper for in-memory queries
- Add memory-mapping layer for persistence
- Add indexing step (sort + build)
- **Pros:** Pure Rust, excellent performance
- **Cons:** Need to implement persistence layer

#### Option C: Port mmmulti::iitree algorithm directly (HIGH effort)
- Implement implicit binary tree in sorted array
- Match C++ interface exactly
- **Pros:** Complete control, no external deps
- **Cons:** ~500 lines to port, needs careful validation

**Recommendation:** Start with **Option A** (FFI wrapper), validate all algorithms work, then consider Option B for long-term.

**Dependencies:**
- `alignments.cpp`: Uses aln_iitree for match storage
- `links.cpp`: Uses node_iitree and path_iitree
- `transclosure.cpp`: Uses q_iitree and s_iitree
- `gfa.cpp`, `vgp.cpp`: Uses node_iitree and path_iitree

### 📋 mmmulti::map and mmmulti::set

**Priority:** MEDIUM (used mainly by transclosure)

**C++ Implementation:**
- Disk-backed sorted arrays
- Succinct bitvectors for key indexing (sdsl::sd_vector)
- IPS4o parallel sorting
- Atomic queue for concurrent writes

**Approach:**

Same three options as iitree. **Recommendation:** Start with FFI wrapper.

**Alternative:** If transclosure is last to migrate, we can rewrite it to use simpler Rust collections (HashMap, BTreeMap) since it's internal to one algorithm.

### ✅ sdsl (Succinct Data Structure Library)

**Priority:** ~~LOW for now~~ **IMPLEMENTED** (seqindex uses these)

**C++ Usage:**
- `sdsl::sd_vector`: Succinct bit vector with rank/select
- `sdsl::csa_wt`: Compressed suffix array (for sequence name index)

**Rust Equivalents (SELECTED):**
- `fm-index 0.3.0`: FM-index with locate support (replaces sdsl::csa_wt)
- `vers-vecs 1.8.1`: Fast rank/select (replaces sdsl::sd_vector)

**Current Status:**
- ✅ seqindex uses fm-index::FMIndexWithLocate for name lookup (O(n log σ) bits)
- ✅ seqindex uses vers_vecs::RsVec for boundaries (O(m log(N/m)) bits)
- **Complexity matches C++ sdsl exactly**
- Still needed for mmmulti structures (will use same approach or FFI)

## Migration Order

### Phase 2A: Enable Core Algorithms ✅ IN PROGRESS

1. ✅ **Step 24:** seqindex basic implementation (COMPLETE)
2. **Step 25:** seqindex FFI bindings for C++
3. **Step 26:** Update seqindex.cpp to use Rust via FFI
4. **Step 27:** Test existing code still works with Rust seqindex

### Phase 2B: Interval Trees

5. **Step 28:** Create FFI wrapper for mmmulti::iitree
   - Opaque IitreeHandle
   - FFI for add(), overlap(), index()
   - Expose to Rust algorithms

6. **Step 29:** (Optional) Start rust-lapper integration research
   - Benchmark vs mmmulti::iitree
   - Design persistence layer
   - Plan migration if beneficial

### Phase 2C: Algorithm Migration

7. **Step 30:** Migrate alignments.cpp to Rust
   - Uses: seqindex, iitree, PAF parsing
   - ~128 lines of core logic
   - Parallel worker pattern

8. **Step 31:** Migrate links.cpp to Rust
   - Uses: seqindex, iitree, sdsl bitvectors
   - ~71 lines
   - May need to FFI wrap more of mmmulti

9. **Step 32:** Migrate compact.cpp to Rust
   - ~85 lines
   - Graph compaction logic

10. **Step 33:** Migrate GFA and VGP output to Rust
    - gfa.cpp: ~207 lines
    - vgp.cpp: ~169 lines
    - Producer/consumer queues, parallel formatting

11. **Step 34:** Migrate transclosure.cpp to Rust (THE BIG ONE)
    - ~643 lines
    - Heavy use of mmmulti structures
    - DisjointSets, atomic_bitvector, atomic_queue
    - May require multiple sub-steps

12. **Step 35:** Migrate main.cpp to Rust
    - CLI parsing (use `clap` crate)
    - Orchestrates all algorithms
    - Final integration

## Design Patterns

### FFI Bridge Pattern

For structures we want to keep in C++ temporarily:

```rust
// Rust side - opaque handle
pub struct IitreeHandle {
    _private: [u8; 0],
}

#[no_mangle]
pub extern "C" fn iitree_create() -> *mut IitreeHandle;

#[no_mangle]
pub extern "C" fn iitree_add(handle: *mut IitreeHandle,
                               start: u64, end: u64, data: u64);

#[no_mangle]
pub extern "C" fn iitree_overlap(handle: *const IitreeHandle,
                                  start: u64, end: u64,
                                  callback: extern "C" fn(...));
```

```cpp
// C++ side - wrapper
class iitree_rust_wrapper {
    IitreeHandle* handle;
public:
    iitree_rust_wrapper() : handle(iitree_create()) {}
    ~iitree_rust_wrapper() { iitree_free(handle); }
    void add(uint64_t start, uint64_t end, uint64_t data) {
        iitree_add(handle, start, end, data);
    }
    // ...
};
```

### Memory-Mapped Collections Pattern

For large disk-backed structures:

```rust
use memmap2::Mmap;

pub struct MmapVec<T> {
    mmap: Option<Mmap>,
    len: usize,
    _phantom: std::marker::PhantomData<T>,
}

impl<T> MmapVec<T> {
    pub fn from_file(path: &str) -> Result<Self, Error> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let len = mmap.len() / std::mem::size_of::<T>();
        Ok(MmapVec { mmap: Some(mmap), len, _phantom: PhantomData })
    }

    pub fn get(&self, index: usize) -> Option<&T> {
        if index < self.len {
            unsafe {
                let ptr = self.mmap.as_ref()?.as_ptr() as *const T;
                Some(&*ptr.add(index))
            }
        } else {
            None
        }
    }
}
```

### Parallel Processing Pattern

Using rayon for parallelism:

```rust
use rayon::prelude::*;

fn process_sequences_parallel(seqindex: &SeqIndex, num_threads: usize) {
    (0..seqindex.n_seqs())
        .into_par_iter()
        .with_max_len(10000)  // Chunk size
        .for_each(|seq_id| {
            // Process sequence
        });
}
```

## Testing Strategy

### 1. Unit Tests
- Each Rust module has comprehensive tests
- Cover edge cases, empty inputs, large inputs
- Current: 93 tests passing

### 2. FFI Integration Tests
- Test C++ can call Rust and vice versa
- Verify memory management (no leaks, no double-frees)
- Test error propagation across FFI boundary

### 3. Algorithmic Validation
- Compare Rust output vs C++ output on same inputs
- Use small test datasets for exact comparison
- Use real datasets for performance benchmarking

### 4. Performance Benchmarks
- Track performance at each migration step
- Identify regressions early
- Use `criterion` crate for Rust benchmarks

## Risk Mitigation

### Risk: Performance Regression
**Mitigation:**
- Benchmark at every step
- Keep C++ version available for comparison
- Profile before optimizing

### Risk: FFI Complexity
**Mitigation:**
- Use cbindgen for automatic header generation
- Minimize FFI surface area
- Clear ownership rules (Rust owns, C++ borrows or vice versa)

### Risk: Memory Leaks Across FFI
**Mitigation:**
- Use valgrind/sanitizers
- Clear ownership documentation
- RAII wrappers on both sides

### Risk: Behavioral Differences
**Mitigation:**
- Extensive test suite
- Diff outputs from C++ and Rust versions
- Gradual rollout (keep both implementations until validated)

## Success Criteria

### Phase 2A Complete (seqindex working):
- ✅ Rust seqindex parses FASTA/FASTQ
- ✅ All query operations implemented
- ⏳ FFI bindings allow C++ to use it
- ⏳ Existing tests pass with Rust seqindex
- ⏳ Performance within 10% of C++ version

### Phase 2B Complete (interval trees working):
- [ ] Rust or FFI-wrapped iitree available
- [ ] All overlap queries work correctly
- [ ] Persistence (save/load) works
- [ ] Performance validated

### Phase 2C Complete (algorithms migrated):
- [ ] All algorithms (alignments, links, compact, transclosure) in Rust
- [ ] GFA/VGP output identical to C++ version
- [ ] End-to-end tests pass
- [ ] Performance meets or exceeds C++ version

### Phase 2 Complete (pure Rust):
- [ ] All C++ code removed (or minimal FFI shims)
- [ ] Single-language codebase
- [ ] All tests passing
- [ ] Documentation updated
- [ ] Ready for Phase 3 (optimization & features)

## Timeline Estimate

- **Phase 2A (seqindex):** 1-2 sessions (mostly done)
- **Phase 2B (interval trees):** 2-3 sessions
- **Phase 2C (algorithms):** 5-8 sessions
- **Total:** 8-13 sessions

## Open Questions

1. ~~**Should we compress sequence names?**~~ **RESOLVED**
   - HashMap uses ~150 bytes/sequence
   - CSA uses ~10 bytes/sequence (O(n log σ) bits)
   - For 10K sequences: 1.5MB vs 100KB
   - **Decision:** ✅ Use CSA (fm-index) - Big O bounds are non-negotiable

2. ~~**Should we use succinct bitvectors now or later?**~~ **RESOLVED**
   - Vec<u64> for offsets uses 8 bytes per sequence
   - sd_vector uses ~0.1 bits per sequence (O(m log(N/m)) bits)
   - For 10K sequences: 80KB vs 1.25KB
   - **Decision:** ✅ Use succinct bitvector (vers-vecs) - Big O bounds are non-negotiable

3. **Pure Rust vs FFI for interval trees?**
   - FFI is faster to implement
   - Pure Rust is better long-term
   - **Decision:** Start FFI, evaluate rust-lapper when algorithms are all ported

4. **Parallel sorting: rayon vs IPS4o?**
   - rayon's par_sort() is good for most cases
   - IPS4o is specialized for huge arrays
   - **Decision:** Try rayon first, benchmark, consider ips4o FFI wrapper if needed

## Next Steps

1. **Immediate:** Add seqindex FFI bindings (Step 25)
2. **Next session:** Update seqindex.cpp to use Rust, test integration (Steps 26-27)
3. **Following session:** Design iitree FFI wrapper (Step 28)

## References

- Phase 2 Research: `PHASE2_RESEARCH.md`
- Rust implementation: `seqwish-rs/src/seqindex.rs` (574 lines, Step 24 complete)
- C++ original: `src/seqindex.cpp`, `src/seqindex.hpp`
- Big O requirements: `.claude/CLAUDE.md`
- Test results: 93/93 passing
