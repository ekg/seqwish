/// Hybrid approach matching C++ behavior:
/// - Use aligned u128 for storage (compiler generates atomic SSE loads)
/// - Use portable_atomic only for CAS operations
/// - This matches the C++ dset64-gccAtomic.hpp approach

use portable_atomic::{AtomicU128, Ordering};

/// 16-byte aligned u128 - compiler will use atomic SSE loads
#[repr(C, align(16))]
struct AlignedU128(u128);

pub struct DisjointSetsAsm {
    data: *mut AlignedU128,
    len: usize,
}

const PARENT_MASK: u128 = u64::MAX as u128;
const RANK_MASK: u128 = (u64::MAX as u128) << 64;

unsafe impl Send for DisjointSetsAsm {}
unsafe impl Sync for DisjointSetsAsm {}

impl DisjointSetsAsm {
    pub fn new(size: usize) -> Self {
        let layout = std::alloc::Layout::array::<AlignedU128>(size).unwrap();
        let ptr = unsafe { std::alloc::alloc(layout) as *mut AlignedU128 };

        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }

        // Initialize
        for i in 0..size {
            unsafe {
                ptr.add(i).write(AlignedU128(i as u128));
            }
        }

        DisjointSetsAsm {
            data: ptr,
            len: size,
        }
    }

    #[inline(always)]
    pub fn find(&self, mut id: usize) -> usize {
        unsafe {
            while id != self.parent_unchecked(id) {
                // Regular load - compiler will use atomic SSE instruction due to alignment
                let value = (*self.data.add(id)).0;
                let new_parent = self.parent_unchecked((value & PARENT_MASK) as usize);
                let new_value = (value & RANK_MASK) | (new_parent as u128);

                if value != new_value {
                    // Use atomic CAS for updates only
                    self.compare_exchange_u128(self.data.add(id) as *mut u128, value, new_value);
                }
                id = new_parent;
            }
            id
        }
    }

    /// Atomic CAS using portable_atomic (matches C++ __sync_bool_compare_and_swap)
    #[inline(always)]
    unsafe fn compare_exchange_u128(&self, ptr: *mut u128, expected: u128, new: u128) -> bool {
        let atomic_ptr = ptr as *const AtomicU128;
        (*atomic_ptr).compare_exchange_weak(
            expected,
            new,
            Ordering::Release,
            Ordering::Acquire,
        ).is_ok()
    }

    #[inline(always)]
    pub fn unite(&self, mut id1: usize, mut id2: usize) -> usize {
        loop {
            id1 = self.find(id1);
            id2 = self.find(id2);

            if id1 == id2 {
                return id1;
            }

            let mut r1 = unsafe { self.rank_unchecked(id1) };
            let mut r2 = unsafe { self.rank_unchecked(id2) };

            if r1 > r2 || (r1 == r2 && id1 < id2) {
                std::mem::swap(&mut r1, &mut r2);
                std::mem::swap(&mut id1, &mut id2);
            }

            let old_entry = ((r1 as u128) << 64) | (id1 as u128);
            let new_entry = ((r1 as u128) << 64) | (id2 as u128);

            unsafe {
                if !self.compare_exchange_u128(self.data.add(id1) as *mut u128, old_entry, new_entry) {
                    continue;
                }

                if r1 == r2 {
                    let old_entry = ((r2 as u128) << 64) | (id2 as u128);
                    let new_entry = (((r2 + 1) as u128) << 64) | (id2 as u128);
                    self.compare_exchange_u128(self.data.add(id2) as *mut u128, old_entry, new_entry);
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
    unsafe fn rank_unchecked(&self, id: usize) -> u64 {
        // Regular load - compiler uses atomic SSE instruction
        let value = (*self.data.add(id)).0;
        ((value >> 64) & PARENT_MASK) as u64
    }

    #[inline(always)]
    unsafe fn parent_unchecked(&self, id: usize) -> usize {
        // Regular load - compiler uses atomic SSE instruction
        let value = (*self.data.add(id)).0;
        (value & PARENT_MASK) as usize
    }
}

impl Drop for DisjointSetsAsm {
    fn drop(&mut self) {
        if !self.data.is_null() {
            unsafe {
                let layout = std::alloc::Layout::array::<AlignedU128>(self.len).unwrap();
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
        let dsets = DisjointSetsAsm::new(10);
        assert_eq!(dsets.find(0), 0);
        dsets.unite(0, 5);
        assert_eq!(dsets.find(0), dsets.find(5));
    }
}
