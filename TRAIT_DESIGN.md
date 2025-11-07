# Trait-Based Interval Tree Architecture

## The Problem

Current implementation is tightly coupled to disk-backed iitree:
- Can't switch to in-memory for small datasets
- Can't benchmark different implementations
- No flexibility for users

## The Solution: Abstract with Traits

```rust
/// Generic interval tree interface
pub trait IntervalTree<K: Ord, V: Clone>: Send + Sync {
    /// Add an interval [start, end) with associated value
    fn add(&mut self, start: K, end: K, value: V) -> io::Result<()>;

    /// Query all intervals that overlap with position
    fn query(&self, pos: K) -> Vec<V>;

    /// Query all intervals that overlap with range [start, end)
    fn query_range(&self, start: K, end: K) -> Vec<V>;

    /// Finalize/index the tree (for disk-backed implementations)
    fn finalize(&mut self) -> io::Result<()>;
}

/// Disk-backed implementation (current iitree-rs)
pub struct DiskBackedTree<K, V> {
    inner: iitree_rs::IITree<K, V>,
}

/// Pure in-memory implementation
pub struct InMemoryTree<K, V> {
    intervals: Vec<(K, K, V)>,  // (start, end, value)
    indexed: bool,
}

/// Smart wrapper that chooses implementation based on size estimate
pub enum AdaptiveTree<K, V> {
    Memory(InMemoryTree<K, V>),
    Disk(DiskBackedTree<K, V>),
}
```

## Usage Pattern

```rust
// At compile time (zero-cost)
fn process_with_tree<T: IntervalTree<u64, PosT>>(
    tree: &mut T,
    data: &[Interval]
) {
    for interval in data {
        tree.add(interval.start, interval.end, interval.value)?;
    }
    tree.finalize()?;

    // Query
    let results = tree.query(position);
}

// At runtime (small overhead, but flexible)
let tree: Box<dyn IntervalTree<u64, PosT>> = if size_estimate < IN_MEMORY_THRESHOLD {
    Box::new(InMemoryTree::new())
} else {
    Box::new(DiskBackedTree::new(path)?)
};

// Or adaptive (best of both)
let tree = AdaptiveTree::with_size_hint(size_estimate);
```

## Implementation Strategy

### Phase 1: Extract Trait
1. Define `IntervalTree` trait
2. Wrap existing iitree-rs in trait impl
3. No behavior change, just abstraction

### Phase 2: Add In-Memory Impl
1. Implement `InMemoryTree` using Vec or BTreeMap
2. Optimize for query performance
3. Compare performance

### Phase 3: Make it Switchable
1. Add command-line flag: `--in-memory`
2. Add size-based auto-detection
3. Add adaptive implementation

## Benefits

1. **Flexibility**: Users can choose based on their dataset
2. **Performance**: Can optimize each implementation separately
3. **Testing**: Easy to benchmark different approaches
4. **Future**: Can add more implementations (e.g., compressed, GPU-backed)

## Potential Performance Wins

### In-Memory Advantages:
- No system calls
- Better cache locality
- Can use SIMD for queries
- Rust Vec optimizations

### Disk-Backed Advantages:
- Handles datasets larger than RAM
- Memory pressure is lower
- Can be persistent between runs

## Key Insight

> The trait abstraction is ZERO COST if used with generics and monomorphization.
> Rust's strength is exactly this: abstraction without overhead.

This is what makes Rust potentially BETTER than C++ - we can have both flexibility
and performance, while C++ would need templates (compile-time only) or virtual
functions (runtime overhead).
