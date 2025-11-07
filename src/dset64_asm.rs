/// ULTRA-UNSAFE version using inline assembly to match C++ __sync builtins
/// This should have ZERO overhead - direct cmpxchg16b inline

use std::arch::asm;

pub struct DisjointSetsAsm {
    data: *mut u128,
    len: usize,
}

const PARENT_MASK: u128 = u64::MAX as u128;
const RANK_MASK: u128 = (u64::MAX as u128) << 64;

unsafe impl Send for DisjointSetsAsm {}
unsafe impl Sync for DisjointSetsAsm {}

impl DisjointSetsAsm {
    pub fn new(size: usize) -> Self {
        let layout = std::alloc::Layout::array::<u128>(size).unwrap();
        let ptr = unsafe { std::alloc::alloc(layout) as *mut u128 };

        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }

        // Initialize
        for i in 0..size {
            unsafe {
                ptr.add(i).write(i as u128);
            }
        }

        DisjointSetsAsm {
            data: ptr,
            len: size,
        }
    }

    /// Compare-exchange using inline assembly (like GCC __sync_bool_compare_and_swap)
    #[inline(always)]
    unsafe fn compare_exchange_u128(
        &self,
        ptr: *mut u128,
        expected: u128,
        new: u128,
    ) -> bool {
        let success: u8;
        let mut expected_lo = expected as u64;
        let mut expected_hi = (expected >> 64) as u64;
        let new_lo = new as u64;
        let new_hi = (new >> 64) as u64;

        // cmpxchg16b requires rbx, so we preserve/restore it
        asm!(
            "xchg rbx, {new_lo}",
            "lock cmpxchg16b [{ptr}]",
            "mov {success}, 0",
            "sete {success}",
            "xchg rbx, {new_lo}",
            ptr = in(reg) ptr,
            new_lo = inout(reg) new_lo => _,
            inout("rax") expected_lo,
            inout("rdx") expected_hi,
            in("rcx") new_hi,
            success = out(reg_byte) success,
        );

        success != 0
    }

    #[inline(always)]
    pub fn find(&self, mut id: usize) -> usize {
        unsafe {
            while id != self.parent_unchecked(id) {
                let ptr = self.data.add(id);
                // Use direct read (not volatile) - cmpxchg provides synchronization
                // This matches C++ behavior and allows compiler optimization
                let value = ptr.read();
                let new_parent = self.parent_unchecked((value & PARENT_MASK) as usize);
                let new_value = (value & RANK_MASK) | (new_parent as u128);

                if value != new_value {
                    // Inline CAS - no function call!
                    self.compare_exchange_u128(ptr, value, new_value);
                }
                id = new_parent;
            }
            id
        }
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
                let ptr = self.data.add(id1);

                // Direct inline assembly - ZERO overhead!
                if !self.compare_exchange_u128(ptr, old_entry, new_entry) {
                    continue;
                }

                if r1 == r2 {
                    let ptr2 = self.data.add(id2);
                    let old_entry = ((r2 as u128) << 64) | (id2 as u128);
                    let new_entry = (((r2 + 1) as u128) << 64) | (id2 as u128);
                    self.compare_exchange_u128(ptr2, old_entry, new_entry);
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
        let ptr = self.data.add(id);
        // Use direct read (not volatile) for better performance
        let value = ptr.read();
        ((value >> 64) & PARENT_MASK) as u64
    }

    #[inline(always)]
    unsafe fn parent_unchecked(&self, id: usize) -> usize {
        let ptr = self.data.add(id);
        // Use direct read (not volatile) for better performance
        let value = ptr.read();
        (value & PARENT_MASK) as usize
    }
}

impl Drop for DisjointSetsAsm {
    fn drop(&mut self) {
        if !self.data.is_null() {
            unsafe {
                let layout = std::alloc::Layout::array::<u128>(self.len).unwrap();
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
