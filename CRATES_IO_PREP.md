# Seqwish Crates.io Preparation

## Checklist for Publishing

### Required ✅

- [x] Update `Cargo.toml` with metadata
  - [x] Change name to `seqwish` (from `seqwish_rs`)
  - [x] Add description, authors, license
  - [x] Add keywords, categories
  - [x] Add repository, readme
- [ ] Create comprehensive README.md
- [ ] Add LICENSE file (MIT)
- [ ] Document all public APIs
- [ ] Add examples directory
- [ ] Add lib.rs documentation
- [ ] Resolve git dependency (iitree-rs needs to be on crates.io first)

### Recommended 📚

- [ ] Add CHANGELOG.md
- [ ] Add more unit tests
- [ ] Add integration test examples
- [ ] Performance benchmarks
- [ ] Usage examples in docs

### Blocking Issues 🚧

1. **iitree-rs dependency**: Currently uses git dependency
   - **Solution**: Either publish iitree-rs first, or inline it, or make it optional

2. **FFI exports**: Currently exports C ABI
   - **Solution**: Make FFI optional feature flag

3. **Binary name**: Currently "seqwish" binary
   - **Already good** - standard practice

## API Documentation Strategy

### Public Library API

The library should expose:

```rust
// Core types
pub use seqindex::SeqIndex;
pub use alignments::unpack_paf_alignments;
pub use transclosure::compute_transitive_closures;
pub use compact::compact_nodes;
pub use links::{derive_links, LinkSet, RankSelectBitVector};
pub use gfa::emit_gfa;

// Configuration
pub struct SeqwishConfig {
    pub num_threads: usize,
    pub repeat_max: u64,
    pub min_repeat_dist: u64,
    pub transclose_batch_size: u64,
    pub show_progress: bool,
}

// High-level API
pub fn build_graph(
    seq_file: &str,
    paf_file: &str,
    config: SeqwishConfig,
) -> Result<VariationGraph, Error>;
```

### Documentation Requirements

1. **Module-level docs** for each major component
2. **Example usage** in lib.rs
3. **API examples** in examples/ directory
4. **Inline examples** in function docs

## RAM-Backed Option Design

### Current Architecture (Disk-Backed)

```
IITree<K, V> → Memory-mapped file on disk
  ├── Scales to huge genomes (100+ GB)
  ├── Slower due to disk I/O
  └── Requires temp directory
```

### Proposed RAM-Backed Option

```rust
// Feature flag approach
[features]
default = ["disk-backed"]
disk-backed = []
ram-backed = []

// Type abstraction
pub trait StorageBackend<K, V> {
    fn new(path: &str) -> io::Result<Self>;
    fn add(&mut self, start: K, end: K, data: V) -> io::Result<()>;
    fn overlap<F>(&self, start: K, end: K, f: F) -> io::Result<()>
    where F: FnMut(K, K, V);
    fn index(&mut self) -> io::Result<()>;
}

// Disk-backed implementation
pub struct DiskIITree<K, V> {
    tree: iitree_rs::IITree<K, V>
}

// RAM-backed implementation
pub struct MemoryIITree<K, V> {
    intervals: Vec<Interval<K, V>>,
    indexed: bool,
}

// Generic over storage
pub struct AlignmentIndex<S: StorageBackend<u64, PosT>> {
    storage: S,
}
```

### Implementation Plan

#### Phase 1: Abstraction Layer
1. Create `StorageBackend` trait
2. Implement for existing `DiskIITree` (thin wrapper)
3. Add generic parameters to all pipeline functions

#### Phase 2: RAM Implementation
1. Implement `MemoryIITree` using `Vec<Interval<K,V>>`
2. Add sorting and binary search for queries
3. Benchmark against disk-backed

#### Phase 3: Configuration
1. Add `--backend [disk|ram]` CLI flag
2. Add `SeqwishConfig::backend` field
3. Use type parameter to select at runtime

#### Phase 4: Optimization
1. Profile RAM-backed version
2. Consider using btree or better structure
3. Add parallel indexing for RAM version

### Trade-offs

**Disk-Backed (Current):**
- ✅ Scales to massive genomes
- ✅ Low memory usage
- ❌ Slower I/O
- ❌ Requires temp space

**RAM-Backed (Proposed):**
- ✅ Faster queries
- ✅ No temp files needed
- ✅ Better for small/medium genomes
- ❌ Limited by RAM
- ❌ Not suitable for huge datasets

### Configuration API

```rust
use seqwish::{SeqwishConfig, StorageBackend};

// Disk-backed (default)
let config = SeqwishConfig::default();

// RAM-backed
let config = SeqwishConfig {
    backend: StorageBackend::Memory,
    ..Default::default()
};

// Auto-select based on size
let config = SeqwishConfig {
    backend: StorageBackend::Auto { threshold_gb: 10 },
    ..Default::default()
};
```

### Estimated Effort

- **Abstraction Layer**: 2-3 days
- **RAM Implementation**: 3-4 days
- **Integration & Testing**: 2-3 days
- **Optimization**: 2-3 days

**Total**: ~2 weeks of focused work

## Next Steps

### Immediate (This Week)
1. ✅ Update CI for Rust-first
2. ✅ Update Cargo.toml metadata
3. [ ] Add comprehensive README
4. [ ] Add LICENSE file
5. [ ] Basic API documentation

### Short-term (Next 2 Weeks)
1. [ ] Work with iitree-rs maintainer to publish it
2. [ ] Add examples directory
3. [ ] Improve test coverage
4. [ ] Add benchmarks

### Medium-term (Next Month)
1. [ ] Implement RAM-backed option
2. [ ] Add feature flags for optional components
3. [ ] Performance optimization
4. [ ] Publish to crates.io

### Long-term
1. [ ] WebAssembly support?
2. [ ] Python bindings?
3. [ ] Incremental graph updates?
4. [ ] Streaming graph construction?
