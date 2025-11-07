//! Generic interval tree abstraction
//!
//! Provides a trait that allows switching between disk-backed and in-memory
//! implementations based on dataset size or user preference.

use std::io;

/// Generic interval tree interface
///
/// Implementations can be disk-backed (for large datasets) or in-memory (for speed).
pub trait IntervalTree<K, V>: Send + Sync
where
    K: Ord + Copy,
    V: Copy,
{
    /// Open writer (for disk-backed trees that use background threads)
    fn open_writer(&mut self) -> io::Result<()> {
        Ok(()) // Default no-op for in-memory trees
    }

    /// Add an interval [start, end) with associated value
    fn add(&mut self, start: K, end: K, value: V) -> io::Result<()>;

    /// Close writer (for disk-backed trees with background threads)
    fn close_writer(&mut self) -> io::Result<()> {
        Ok(()) // Default no-op for in-memory trees
    }

    /// Query all intervals that overlap with range [start, end)
    /// Calls the provided function for each overlapping interval with (index, start, end, value)
    fn overlap<F>(&self, start: K, end: K, func: F) -> io::Result<()>
    where
        F: FnMut(usize, K, K, V);

    /// Finalize the tree (sort, index, etc.) before querying
    /// Also available as `index()` for compatibility
    fn finalize(&mut self) -> io::Result<()>;

    /// Index the tree (alias for finalize, for IITree compatibility)
    fn index(&mut self) -> io::Result<()> {
        self.finalize()
    }

    /// Get the number of intervals stored
    fn len(&self) -> usize;

    /// Check if the tree is empty
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Disk-backed interval tree using iitree-rs
///
/// Stores all data on disk via memory mapping. Good for datasets that
/// don't fit in RAM or when you want persistent storage.
pub mod disk {
    use super::*;
    use iitree_rs::IITree;
    use std::path::Path;

    pub struct DiskBackedTree<K, V>
    where
        K: Ord + Copy + Default + Send + Sync + 'static,
        V: Copy + Default + Send + Sync + 'static,
    {
        inner: IITree<K, V>,
    }

    impl<K, V> DiskBackedTree<K, V>
    where
        K: Ord + Copy + Default + Send + Sync + 'static,
        V: Copy + Default + Send + Sync + 'static,
    {
        pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
            Ok(DiskBackedTree {
                inner: IITree::new(path)?,
            })
        }
    }

    impl<K, V> IntervalTree<K, V> for DiskBackedTree<K, V>
    where
        K: Ord + Copy + Default + Send + Sync + 'static,
        V: Copy + Default + Send + Sync + 'static,
    {
        fn open_writer(&mut self) -> io::Result<()> {
            self.inner.open_writer()
        }

        fn add(&mut self, start: K, end: K, value: V) -> io::Result<()> {
            self.inner.add(start, end, value)
        }

        fn close_writer(&mut self) -> io::Result<()> {
            self.inner.close_writer()
        }

        fn overlap<F>(&self, start: K, end: K, func: F) -> io::Result<()>
        where
            F: FnMut(usize, K, K, V),
        {
            self.inner.overlap(start, end, func)
        }

        fn finalize(&mut self) -> io::Result<()> {
            self.inner.index()
        }

        fn len(&self) -> usize {
            self.inner.len()
        }
    }
}

/// In-memory interval tree using Vec
///
/// Stores all intervals in RAM. Much faster than disk-backed for datasets
/// that fit in memory. Uses simple binary search after sorting.
pub mod memory {
    use super::*;
    use rayon::prelude::*;

    #[derive(Clone)]
    struct Interval<K, V> {
        start: K,
        end: K,
        max: K, // Maximum end position in subtree (for augmented tree)
        value: V,
    }

    pub struct InMemoryTree<K, V>
    where
        K: Ord + Copy,
        V: Copy,
    {
        intervals: Vec<Interval<K, V>>,
        indexed: bool,
        max_level: usize, // Maximum level in the implicit binary tree
    }

    impl<K, V> InMemoryTree<K, V>
    where
        K: Ord + Copy,
        V: Copy,
    {
        pub fn new() -> Self {
            InMemoryTree {
                intervals: Vec::new(),
                indexed: false,
                max_level: 0,
            }
        }

        pub fn with_capacity(capacity: usize) -> Self {
            InMemoryTree {
                intervals: Vec::with_capacity(capacity),
                indexed: false,
                max_level: 0,
            }
        }

        /// IAITree overlap query with stack-based traversal
        /// This is a direct port of the iitree-rs overlap algorithm
        /// Uses the augmented max values to prune entire subtrees
        fn overlap_iaitree<F>(&self, query_start: K, query_end: K, mut func: F)
        where
            F: FnMut(usize, K, K, V),
        {
            let n = self.intervals.len();
            if n == 0 {
                return;
            }

            #[derive(Clone, Copy)]
            struct StackCell {
                k: usize, // level
                x: usize, // node index
                w: u8,    // 0 if left child not processed, 1 if processed
            }

            let mut stack = [StackCell { k: 0, x: 0, w: 0 }; 64];
            let mut t = 0; // stack pointer

            // Push root
            stack[t] = StackCell {
                k: self.max_level,
                x: (1 << self.max_level) - 1,
                w: 0,
            };
            t += 1;

            // Top-down traversal
            while t > 0 {
                t -= 1;
                let z = stack[t];

                if z.k <= 3 {
                    // Small subtree: linear scan
                    let i0 = (z.x >> z.k) << z.k;
                    let i1 = std::cmp::min(i0 + (1 << (z.k + 1)) - 1, n);

                    for i in i0..i1 {
                        if self.intervals[i].start >= query_end {
                            break;
                        }
                        if query_start < self.intervals[i].end {
                            func(
                                i,
                                self.intervals[i].start,
                                self.intervals[i].end,
                                self.intervals[i].value,
                            );
                        }
                    }
                } else if z.w == 0 {
                    // Left child not processed
                    let y = z.x - (1 << (z.k - 1));

                    // Re-add this node with left child marked as processed
                    stack[t] = StackCell {
                        k: z.k,
                        x: z.x,
                        w: 1,
                    };
                    t += 1;

                    // Push left child if it might overlap
                    if y >= n || self.intervals[y].max > query_start {
                        stack[t] = StackCell {
                            k: z.k - 1,
                            x: y,
                            w: 0,
                        };
                        t += 1;
                    }
                } else if z.x < n && self.intervals[z.x].start < query_end {
                    // Process current node and push right child
                    if query_start < self.intervals[z.x].end {
                        func(
                            z.x,
                            self.intervals[z.x].start,
                            self.intervals[z.x].end,
                            self.intervals[z.x].value,
                        );
                    }

                    // Push right child
                    stack[t] = StackCell {
                        k: z.k - 1,
                        x: z.x + (1 << (z.k - 1)),
                        w: 0,
                    };
                    t += 1;
                }
            }
        }
    }

    impl<K, V> IntervalTree<K, V> for InMemoryTree<K, V>
    where
        K: Ord + Copy + Send + Sync,
        V: Copy + Send + Sync,
    {
        fn add(&mut self, start: K, end: K, value: V) -> io::Result<()> {
            self.intervals.push(Interval {
                start,
                end,
                max: end, // Initially set max to end
                value,
            });
            self.indexed = false;
            Ok(())
        }

        fn overlap<F>(&self, query_start: K, query_end: K, func: F) -> io::Result<()>
        where
            F: FnMut(usize, K, K, V),
        {
            if !self.indexed {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Tree not indexed - call finalize() first",
                ));
            }

            self.overlap_iaitree(query_start, query_end, func);
            Ok(())
        }

        fn finalize(&mut self) -> io::Result<()> {
            if self.indexed {
                return Ok(());
            }

            let n = self.intervals.len();
            if n == 0 {
                self.indexed = true;
                self.max_level = 0;
                return Ok(());
            }

            // Sort by start position, then by end (matching iitree-rs behavior)
            if n > 10000 {
                self.intervals.par_sort_unstable_by(|a, b| {
                    a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end))
                });
            } else {
                self.intervals
                    .sort_unstable_by(|a, b| a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end)));
            }

            // Build augmented tree (compute max values bottom-up)
            // This is a direct port of iitree-rs index_core()

            // Initialize leaves (level 0) - every other node starting from 0
            let mut last_i = 0;
            let mut last = self.intervals[0].end;
            for i in (0..n).step_by(2) {
                last_i = i;
                self.intervals[i].max = self.intervals[i].end;
                last = self.intervals[i].max;
            }

            // Process internal nodes bottom-up
            let mut k = 1;
            while (1usize << k) <= n {
                let x = 1usize << (k - 1);
                let i0 = (x << 1) - 1; // First node at this level
                let step = x << 2;

                // Traverse all nodes at level k
                let mut i = i0;
                while i < n {
                    // max value of left child
                    let el = self.intervals[i - x].max;

                    // max value of right child (may not exist)
                    let er = if i + x < n {
                        self.intervals[i + x].max
                    } else {
                        last
                    };

                    // Compute max for this node
                    let mut e = self.intervals[i].end;
                    if el > e {
                        e = el;
                    }
                    if er > e {
                        e = er;
                    }
                    self.intervals[i].max = e;

                    i += step;
                }

                // Update last_i to point to parent of original last_i
                last_i = if (last_i >> k) & 1 != 0 {
                    last_i - x
                } else {
                    last_i + x
                };
                if last_i < n {
                    last = self.intervals[last_i].max;
                }

                k += 1;
            }

            self.max_level = k - 1;
            self.indexed = true;
            Ok(())
        }

        fn len(&self) -> usize {
            self.intervals.len()
        }
    }

    impl<K, V> Default for InMemoryTree<K, V>
    where
        K: Ord + Copy,
        V: Copy,
    {
        fn default() -> Self {
            Self::new()
        }
    }
}

/// Adaptive tree that can be either disk-backed or in-memory
///
/// This enum allows runtime selection of implementation while maintaining
/// good performance through match-based dispatch (no virtual function overhead).
pub enum AdaptiveTree<K, V>
where
    K: Ord + Copy + Default + Send + Sync + 'static,
    V: Copy + Default + Send + Sync + 'static,
{
    Disk(disk::DiskBackedTree<K, V>),
    Memory(memory::InMemoryTree<K, V>),
}

impl<K, V> AdaptiveTree<K, V>
where
    K: Ord + Copy + Default + Send + Sync + 'static,
    V: Copy + Default + Send + Sync + 'static,
{
    /// Create a disk-backed tree
    pub fn new_disk<P: AsRef<std::path::Path>>(path: P) -> io::Result<Self> {
        Ok(AdaptiveTree::Disk(disk::DiskBackedTree::new(path)?))
    }

    /// Create a true in-memory tree (no filesystem dependency)
    pub fn new_memory() -> io::Result<Self> {
        Ok(AdaptiveTree::Memory(memory::InMemoryTree::new()))
    }

    /// Create an in-memory tree with capacity hint
    pub fn new_memory_with_capacity(capacity: usize) -> io::Result<Self> {
        Ok(AdaptiveTree::Memory(memory::InMemoryTree::with_capacity(
            capacity,
        )))
    }
}

impl<K, V> IntervalTree<K, V> for AdaptiveTree<K, V>
where
    K: Ord + Copy + Default + Send + Sync + 'static,
    V: Copy + Default + Send + Sync + 'static,
{
    fn open_writer(&mut self) -> io::Result<()> {
        match self {
            AdaptiveTree::Disk(tree) => tree.open_writer(),
            AdaptiveTree::Memory(tree) => tree.open_writer(),
        }
    }

    fn add(&mut self, start: K, end: K, value: V) -> io::Result<()> {
        match self {
            AdaptiveTree::Disk(tree) => tree.add(start, end, value),
            AdaptiveTree::Memory(tree) => tree.add(start, end, value),
        }
    }

    fn close_writer(&mut self) -> io::Result<()> {
        match self {
            AdaptiveTree::Disk(tree) => tree.close_writer(),
            AdaptiveTree::Memory(tree) => tree.close_writer(),
        }
    }

    fn overlap<F>(&self, start: K, end: K, func: F) -> io::Result<()>
    where
        F: FnMut(usize, K, K, V),
    {
        match self {
            AdaptiveTree::Disk(tree) => tree.overlap(start, end, func),
            AdaptiveTree::Memory(tree) => tree.overlap(start, end, func),
        }
    }

    fn finalize(&mut self) -> io::Result<()> {
        match self {
            AdaptiveTree::Disk(tree) => tree.finalize(),
            AdaptiveTree::Memory(tree) => tree.finalize(),
        }
    }

    fn len(&self) -> usize {
        match self {
            AdaptiveTree::Disk(tree) => tree.len(),
            AdaptiveTree::Memory(tree) => tree.len(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tempfile;

    #[test]
    fn test_in_memory_tree() {
        let mut tree = memory::InMemoryTree::new();

        tree.add(10, 20, 1).unwrap();
        tree.add(15, 25, 2).unwrap();
        tree.add(30, 40, 3).unwrap();

        tree.finalize().unwrap();

        // Query range [17, 18) should find intervals 1 and 2
        let mut results = Vec::new();
        tree.overlap(17, 18, |_idx, _start, _end, value| {
            results.push(value);
        })
        .unwrap();
        assert_eq!(results.len(), 2);
        assert!(results.contains(&1));
        assert!(results.contains(&2));

        // Query range [35, 36) should find interval 3
        let mut results = Vec::new();
        tree.overlap(35, 36, |_idx, _start, _end, value| {
            results.push(value);
        })
        .unwrap();
        assert_eq!(results, vec![3]);

        // Query range [5, 6) should find nothing
        let mut results = Vec::new();
        tree.overlap(5, 6, |_idx, _start, _end, value| {
            results.push(value);
        })
        .unwrap();
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_memory_vs_disk_equivalence() {
        use std::collections::HashSet;

        // Create both tree types
        let mut mem_tree = memory::InMemoryTree::new();
        let disk_path = tempfile::create("test-equiv-", ".iit").unwrap();
        let mut disk_tree = disk::DiskBackedTree::new(&disk_path).unwrap();

        // Add same intervals to both - including ones with same start but different ends
        let intervals = vec![
            (10, 20, 1u64),
            (10, 25, 2u64), // Same start, different end
            (10, 15, 3u64), // Same start, different end
            (15, 25, 4u64),
            (30, 40, 5u64),
        ];

        disk_tree.open_writer().unwrap();
        for (start, end, val) in &intervals {
            mem_tree.add(*start, *end, *val).unwrap();
            disk_tree.add(*start, *end, *val).unwrap();
        }
        disk_tree.close_writer().unwrap();

        mem_tree.finalize().unwrap();
        disk_tree.finalize().unwrap();

        // Query both and collect results
        let mut mem_results = Vec::new();
        let mut disk_results = Vec::new();

        mem_tree
            .overlap(12, 18, |idx, start, end, val| {
                mem_results.push((idx, start, end, val));
            })
            .unwrap();

        disk_tree
            .overlap(12, 18, |idx, start, end, val| {
                disk_results.push((idx, start, end, val));
            })
            .unwrap();

        // Compare ignoring index (since that might differ due to sort stability)
        let mem_set: HashSet<_> = mem_results.iter().map(|(_, s, e, v)| (s, e, v)).collect();
        let disk_set: HashSet<_> = disk_results.iter().map(|(_, s, e, v)| (s, e, v)).collect();

        assert_eq!(
            mem_set, disk_set,
            "Memory and disk trees should return the same intervals.\n\
             Memory: {:?}\n\
             Disk: {:?}",
            mem_results, disk_results
        );

        // Clean up temp file
        tempfile::remove(&disk_path);
    }
}
