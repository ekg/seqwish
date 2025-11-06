/// UNSAFE version of dset64 that eliminates bounds checking
/// for maximum performance matching C++
///
/// This matches the C++ implementation which uses raw pointers.

use portable_atomic::{AtomicU128, Ordering};

/// Lock-free disjoint set using raw pointer for zero-overhead array access
pub struct DisjointSetsUnsafe {
    data: *mut AtomicU128,
    len: usize,
}

const PARENT_MASK: u128 = u64::MAX as u128;
const RANK_MASK: u128 = (u64::MAX as u128) << 64;

unsafe impl Send for DisjointSetsUnsafe {}
unsafe impl Sync for DisjointSetsUnsafe {}

impl DisjointSetsUnsafe {
    /// Create a new disjoint set with `size` elements
    pub fn new(size: usize) -> Self {
        assert!(
            AtomicU128::is_lock_free(),
            "AtomicU128 must be lock-free!"
        );

        // Allocate aligned memory like C++
        let layout = std::alloc::Layout::array::<AtomicU128>(size).unwrap();
        let ptr = unsafe { std::alloc::alloc(layout) as *mut AtomicU128 };

        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }

        // Initialize each element to point to itself (parent = self)
        for i in 0..size {
            unsafe {
                ptr.add(i).write(AtomicU128::new(i as u128));
            }
        }

        DisjointSetsUnsafe {
            data: ptr,
            len: size,
        }
    }

    /// Find with NO bounds checking (matches C++ performance)
    #[inline(always)]
    pub fn find(&self, mut id: usize) -> usize {
        unsafe {
            while id != self.parent_unchecked(id) {
                let ptr = self.data.add(id);
                let value = (*ptr).load(Ordering::Relaxed);
                let new_parent = self.parent_unchecked((value & PARENT_MASK) as usize);
                let new_value = (value & RANK_MASK) | (new_parent as u128);

                if value != new_value {
                    let _ = (*ptr).compare_exchange_weak(
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
    }

    /// Unite with NO bounds checking (matches C++ performance)
    #[inline(always)]
    pub fn unite(&self, mut id1: usize, mut id2: usize) -> usize {
        loop {
            id1 = self.find(id1);
            id2 = self.find(id2);

            if id1 == id2 {
                return id1;
            }

            let mut r1 = self.rank_unchecked(id1);
            let mut r2 = self.rank_unchecked(id2);

            if r1 > r2 || (r1 == r2 && id1 < id2) {
                std::mem::swap(&mut r1, &mut r2);
                std::mem::swap(&mut id1, &mut id2);
            }

            let old_entry = ((r1 as u128) << 64) | (id1 as u128);
            let new_entry = ((r1 as u128) << 64) | (id2 as u128);

            unsafe {
                let ptr = self.data.add(id1);
                if (*ptr).compare_exchange(
                    old_entry,
                    new_entry,
                    Ordering::SeqCst,
                    Ordering::SeqCst,
                ).is_err() {
                    continue;
                }

                if r1 == r2 {
                    let ptr2 = self.data.add(id2);
                    let old_entry = ((r2 as u128) << 64) | (id2 as u128);
                    let new_entry = (((r2 + 1) as u128) << 64) | (id2 as u128);
                    let _ = (*ptr2).compare_exchange_weak(
                        old_entry,
                        new_entry,
                        Ordering::SeqCst,
                        Ordering::SeqCst,
                    );
                }

                return id2;
            }
        }
    }

    #[inline(always)]
    pub fn size(&self) -> usize {
        self.len
    }

    #[inline(always)]
    fn rank_unchecked(&self, id: usize) -> u64 {
        unsafe {
            let ptr = self.data.add(id);
            (((*ptr).load(Ordering::Relaxed) >> 64) & PARENT_MASK) as u64
        }
    }

    #[inline(always)]
    fn parent_unchecked(&self, id: usize) -> usize {
        unsafe {
            let ptr = self.data.add(id);
            ((*ptr).load(Ordering::Relaxed) & PARENT_MASK) as usize
        }
    }
}

impl Drop for DisjointSetsUnsafe {
    fn drop(&mut self) {
        if !self.data.is_null() {
            unsafe {
                // Drop each AtomicU128
                for i in 0..self.len {
                    self.data.add(i).drop_in_place();
                }
                // Deallocate
                let layout = std::alloc::Layout::array::<AtomicU128>(self.len).unwrap();
                std::alloc::dealloc(self.data as *mut u8, layout);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let dsets = DisjointSetsUnsafe::new(10);
        assert_eq!(dsets.find(0), 0);
        assert_eq!(dsets.find(5), 5);

        dsets.unite(0, 5);
        assert_eq!(dsets.find(0), dsets.find(5));
    }
}
