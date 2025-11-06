/// Lock-free parallel disjoint set data structure (aka UNION-FIND)
/// with path compression and union by rank
///
/// Supports concurrent find(), same() and unite() calls as described
/// in the paper "Wait-free Parallel Algorithms for the Union-Find Problem"
/// by Richard J. Anderson and Heather Woll
///
/// This is a Rust port of dset64-gccAtomic.hpp from seqwish-cpp.
/// Uses 128-bit atomic primitives for 64-bit item ids.

use portable_atomic::{AtomicU128, Ordering};

/// Lock-free disjoint set data structure
pub struct DisjointSets {
    data: Vec<AtomicU128>,
}

// Masks for extracting parent and rank from 128-bit value
// Lower 64 bits = parent, upper 64 bits = rank
const PARENT_MASK: u128 = u64::MAX as u128;
const RANK_MASK: u128 = (u64::MAX as u128) << 64;

impl DisjointSets {
    /// Create a new disjoint set with `size` elements
    /// Each element starts in its own set (parent = self)
    pub fn new(size: usize) -> Self {
        // Verify that AtomicU128 is lock-free (uses cmpxchg16b on x86_64)
        assert!(
            AtomicU128::is_lock_free(),
            "AtomicU128 must be lock-free for dset64 performance!"
        );

        let data: Vec<AtomicU128> = (0..size)
            .map(|i| AtomicU128::new(i as u128))
            .collect();

        DisjointSets { data }
    }

    /// Find the representative of the set containing `id`
    /// Uses path compression for efficiency
    #[inline(always)]
    pub fn find(&self, mut id: usize) -> usize {
        while id != self.parent(id) {
            let value = self.data[id].load(Ordering::Relaxed);
            let new_parent = self.parent((value & PARENT_MASK) as usize);
            let new_value = (value & RANK_MASK) | (new_parent as u128);

            // Try to update parent (may fail due to concurrent access, that's ok)
            // Use SeqCst to match C++ __sync_bool_compare_and_swap
            if value != new_value {
                let _ = self.data[id].compare_exchange_weak(
                    value,
                    new_value,
                    Ordering::SeqCst,
                    Ordering::SeqCst,
                );
            }
            id = new_parent;
        }
        id
    }

    /// Check if two elements are in the same set
    #[inline]
    pub fn same(&self, mut id1: usize, mut id2: usize) -> bool {
        loop {
            id1 = self.find(id1);
            id2 = self.find(id2);

            if id1 == id2 {
                return true;
            }
            if self.parent(id1) == id1 {
                return false;
            }
        }
    }

    /// Unite the sets containing `id1` and `id2`
    /// Returns the representative of the merged set
    #[inline(always)]
    pub fn unite(&self, mut id1: usize, mut id2: usize) -> usize {
        loop {
            id1 = self.find(id1);
            id2 = self.find(id2);

            if id1 == id2 {
                return id1;
            }

            let mut r1 = self.rank(id1);
            let mut r2 = self.rank(id2);

            // Union by rank: attach smaller rank tree under root of higher rank tree
            if r1 > r2 || (r1 == r2 && id1 < id2) {
                std::mem::swap(&mut r1, &mut r2);
                std::mem::swap(&mut id1, &mut id2);
            }

            let old_entry = ((r1 as u128) << 64) | (id1 as u128);
            let new_entry = ((r1 as u128) << 64) | (id2 as u128);

            // Try to make id1 point to id2
            // Use SeqCst to match C++ __sync_bool_compare_and_swap
            if self.data[id1].compare_exchange(
                old_entry,
                new_entry,
                Ordering::SeqCst,
                Ordering::SeqCst,
            ).is_err() {
                continue; // CAS failed, retry
            }

            // If ranks were equal, increment rank of new root
            if r1 == r2 {
                let old_entry = ((r2 as u128) << 64) | (id2 as u128);
                let new_entry = (((r2 + 1) as u128) << 64) | (id2 as u128);
                // Try to update the rank (may fail, that's ok)
                let _ = self.data[id2].compare_exchange_weak(
                    old_entry,
                    new_entry,
                    Ordering::SeqCst,
                    Ordering::SeqCst,
                );
            }

            return id2;
        }
    }

    /// Get the number of elements
    #[inline]
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// Get the rank of an element
    #[inline(always)]
    fn rank(&self, id: usize) -> u64 {
        ((self.data[id].load(Ordering::Relaxed) >> 64) & PARENT_MASK) as u64
    }

    /// Get the parent of an element
    #[inline(always)]
    fn parent(&self, id: usize) -> usize {
        (self.data[id].load(Ordering::Relaxed) & PARENT_MASK) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let dsets = DisjointSets::new(10);

        // Initially, each element is in its own set
        assert_eq!(dsets.find(0), 0);
        assert_eq!(dsets.find(5), 5);
        assert!(!dsets.same(0, 5));

        // Unite 0 and 5
        dsets.unite(0, 5);
        assert!(dsets.same(0, 5));

        // Unite 5 and 7
        dsets.unite(5, 7);
        assert!(dsets.same(0, 7));
        assert!(dsets.same(5, 7));

        // 3 is still separate
        assert!(!dsets.same(0, 3));
    }

    #[test]
    fn test_path_compression() {
        let dsets = DisjointSets::new(5);

        // Create a chain: 0 -> 1 -> 2 -> 3 -> 4
        dsets.unite(0, 1);
        dsets.unite(1, 2);
        dsets.unite(2, 3);
        dsets.unite(3, 4);

        // All should be in the same set
        let root = dsets.find(0);
        assert_eq!(dsets.find(1), root);
        assert_eq!(dsets.find(2), root);
        assert_eq!(dsets.find(3), root);
        assert_eq!(dsets.find(4), root);
    }
}
