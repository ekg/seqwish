# Transclosure Migration Plan

## Overview
Migrate `transclosure.cpp` (643 lines) to Rust. This module computes transitive closures for variation graph construction.

## Status: 🔴 Not Started

---

## Dependencies Matrix

| C++ Dependency | Rust Equivalent | Status | Notes |
|----------------|-----------------|--------|-------|
| atomic_bitvector | `bitvec` (atomic) | 🟡 TODO | Enable atomic feature |
| atomic_queue | `crossbeam-queue` | 🟡 TODO | Use ArrayQueue |
| dset64-gccAtomic | `uf_rush` | 🟡 TODO | Lock-free union-find |
| spinlock | `parking_lot` | 🟡 TODO | Prefer Mutex over spinlock |
| wang hash | Direct port | 🟡 TODO | Port from wang.hpp |
| **ips4o** | **rayon** | 🟡 TODO | **Use par_sort_unstable** |
| SDSL bit_vector | `simple-sds` | 🟡 TODO | For rank/select operations |
| flat_hash_map | `hashbrown` | ✅ READY | Already in std |
| paryfor | `rayon` | ✅ READY | Use par_iter |

---

## Cargo Dependencies

```toml
[dependencies]
rayon = "1.11"                  # Parallel iteration & sorting
crossbeam-queue = "0.3"         # Lock-free queues
bitvec = { version = "1.0", features = ["atomic"] }
parking_lot = "0.12"            # Mutexes
simple-sds = "0.8"              # Succinct data structures
uf_rush = "0.1"                 # Union-find
```

---

## Migration Checklist

### Phase 1: Setup Dependencies ⏳
- [ ] Add dependencies to `seqwish-rs/Cargo.toml`
- [ ] Port wang_hash_64 function
- [ ] Create type aliases for queues
- [ ] Test dependency compatibility

### Phase 2: Core Functions 🔴
- [ ] `extend_range()` - Extend position ranges
- [ ] `flush_range()` - Write range to iitree
- [ ] `flush_ranges()` - Flush all buffered ranges
- [ ] `for_each_fresh_range()` - Iterate unseen ranges
- [ ] `handle_range()` - Process range overlaps
- [ ] `explore_overlaps()` - Find transitive matches

### Phase 3: Main Algorithm 🔴
- [ ] `write_graph_chunk()` - Write graph sequence
- [ ] `compute_transitive_closures()` - Main entry point
- [ ] Parallel processing with rayon
- [ ] Atomic bitvector operations

### Phase 4: Testing & Validation 🔴
- [ ] Unit tests for helper functions
- [ ] Integration test with small dataset
- [ ] Verify output matches C++ byte-for-byte
- [ ] Performance benchmarks
- [ ] All 30 HLA tests pass

### Phase 5: Integration 🔴
- [ ] Update main.cpp to call Rust version
- [ ] FFI wrapper in lib.rs
- [ ] Remove C++ transclosure code
- [ ] Update build system
- [ ] Commit & push

---

## Key Algorithm Notes

**Transitive Closure**:
- Computes equivalence classes of aligned positions
- Uses interval trees for range queries
- Parallel exploration with work queue
- Union-find for connected components

**Sorting Requirements**:
- Uses parallel sorting on large vectors
- C++ uses ips4o (sample-based, ~3x faster)
- Rust uses rayon pdqsort (quicksort-based, still very fast)
- Performance difference acceptable for migration

**Performance Critical**:
- Atomic bitvector operations (mark visited)
- Queue push/pop in parallel workers
- Union-find unite operations
- Large parallel sorts

---

## Success Criteria

✅ All 30 tests pass with identical MD5 hashes
✅ No performance regression >10%
✅ Pure Rust implementation (no C++ dependencies)
✅ Thread-safe by construction (Rust ownership)

---

## Notes

- Sorting: Decided to use `rayon::par_sort_unstable()` instead of incomplete ips4o port
- Union-find: Using `uf_rush` crate; may need custom impl if performance issues
- SDSL: Using `simple-sds` for rank/select; may need multiple crates for full functionality
- Bitvector: `bitvec` with atomic feature provides thread-safe operations

---

## Estimated Effort

- Setup: 1 hour
- Core functions: 4-6 hours
- Main algorithm: 4-6 hours
- Testing: 2-4 hours
- Integration: 1-2 hours

**Total: 12-19 hours** (1.5-2.5 days)
