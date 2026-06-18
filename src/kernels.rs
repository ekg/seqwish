use cuda_device::atomic::{AtomicOrdering, DeviceAtomicU32};
use cuda_device::{cuda_module, kernel, thread, DisjointSlice};

#[cuda_module]
pub mod kernels {
    use super::*;

    #[kernel]
    pub fn initialize_parents(mut parents: DisjointSlice<u32>, n: u32) {
        // Sets parents[i] = i for all i < n.
        let idx_usize = thread::index_1d().get();
        if idx_usize < n as usize {
            unsafe {
                *parents.as_mut_ptr().add(idx_usize) = idx_usize as u32;
            }
        }
    }

    #[kernel]
    pub fn union_step(mut parents: DisjointSlice<u32>, edges: &[[u32; 2]], num_edges: i32) {
        // Shiloach-Vishkin union step with path-halving.
        //
        // Each unsafe block below is sound because:
        //   - root_u, p_u, gp_u, root_v, p_v, gp_v are always
        //     values read back from the parents array, so they are in range of n.
        //   - The slice was allocated with exactly n elements.
        let edge_idx = thread::index_1d().get();
        if edge_idx >= num_edges as usize {
            return;
        }

        let edge = edges[edge_idx];
        let mut u = edge[0];
        let mut v = edge[1];

        let base = parents.as_mut_ptr();

        loop {
            // --- Find root of u with path-halving ---
            let mut root_u = u;
            loop {
                // SAFETY: root_u was read from parents, so it is a valid index.
                let p_u = unsafe { *base.add(root_u as usize) };
                if p_u == root_u {
                    break;
                }
                let gp_u = unsafe { *base.add(p_u as usize) };
                // Path-halving write, so skip to grandparent.
                unsafe {
                    *base.add(root_u as usize) = gp_u;
                }
                root_u = p_u;
            }

            // --- Find root of v with path-halving ---
            let mut root_v = v;
            loop {
                // SAFETY: same reasoning as root_u.
                let p_v = unsafe { *base.add(root_v as usize) };
                if p_v == root_v {
                    break;
                }
                let gp_v = unsafe { *base.add(p_v as usize) };
                unsafe {
                    *base.add(root_v as usize) = gp_v;
                }
                root_v = p_v;
            }

            if root_u == root_v {
                break;
            }

            // Canonical ordering: smaller root points to larger.
            if root_u > root_v {
                let tmp = root_u;
                root_u = root_v;
                root_v = tmp;
            }

            // Atomic CAS to link root_u -> root_v.
            // SAFETY: root_u is a valid index; DeviceAtomicU32 is repr(transparent).
            let ptr = unsafe { &*(base.add(root_u as usize) as *const DeviceAtomicU32) };
            match ptr.compare_exchange(
                root_u,
                root_v,
                AtomicOrdering::Relaxed,
                AtomicOrdering::Relaxed,
            ) {
                Ok(_) => break,
                Err(old_val) => {
                    u = old_val;
                    v = root_v;
                }
            }
        }
    }

    #[kernel]
    pub fn pointer_jump(mut parents: DisjointSlice<u32>, n: u32, mut changed: DisjointSlice<u32>) {
        // Pointer-jump to flatten the forest in one pass.
        //
        // SAFETY for each block: idx_usize < n (guarded above), and `p` is a
        // value read from the parents array so it is also in [0, n).
        let idx_usize = thread::index_1d().get();
        if idx_usize >= n as usize {
            return;
        }

        let base = parents.as_mut_ptr();
        let p = unsafe { *base.add(idx_usize) };
        let gp = unsafe { *base.add(p as usize) };

        if p != gp {
            unsafe {
                *base.add(idx_usize) = gp;
            }
            // SAFETY: changed has exactly one element.
            let ptr = unsafe { &*(changed.as_mut_ptr() as *const DeviceAtomicU32) };
            ptr.fetch_or(1, AtomicOrdering::Relaxed);
        }
    }
}
