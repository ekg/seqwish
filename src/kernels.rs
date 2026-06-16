use cuda_device::atomic::{AtomicOrdering, DeviceAtomicU32};
use cuda_device::{cuda_module, kernel, thread, DisjointSlice};

#[cuda_module]
pub mod kernels {
    use super::*;

    #[kernel]
    pub fn initialize_parents(mut parents: DisjointSlice<u32>, n: u32) {
        let idx = thread::index_1d();
        let idx_usize = idx.get();
        if idx_usize < n as usize {
            unsafe {
                *parents.as_mut_ptr().add(idx_usize) = idx_usize as u32;
            }
        }
    }

    #[kernel]
    pub fn union_step(mut parents: DisjointSlice<u32>, edges: &[[u32; 2]], num_edges: i32) {
        let edge_idx = thread::index_1d();
        let edge_idx_usize = edge_idx.get();
        if edge_idx_usize >= num_edges as usize {
            return;
        }

        let edge = edges[edge_idx_usize];
        let mut u = edge[0];
        let mut v = edge[1];

        loop {
            let mut root_u = u;
            loop {
                // Safe read using parents pointer
                let p_u = unsafe { *parents.as_mut_ptr().add(root_u as usize) };
                if p_u == root_u {
                    break;
                }
                let gp_u = unsafe { *parents.as_mut_ptr().add(p_u as usize) };
                unsafe {
                    *parents.as_mut_ptr().add(root_u as usize) = gp_u;
                }
                root_u = p_u;
            }

            let mut root_v = v;
            loop {
                let p_v = unsafe { *parents.as_mut_ptr().add(root_v as usize) };
                if p_v == root_v {
                    break;
                }
                let gp_v = unsafe { *parents.as_mut_ptr().add(p_v as usize) };
                unsafe {
                    *parents.as_mut_ptr().add(root_v as usize) = gp_v;
                }
                root_v = p_v;
            }

            if root_u == root_v {
                break;
            }

            if root_u > root_v {
                let tmp = root_u;
                root_u = root_v;
                root_v = tmp;
            }

            // Perform atomic CAS using DeviceAtomicU32
            let ptr =
                unsafe { &*(parents.as_mut_ptr().add(root_u as usize) as *const DeviceAtomicU32) };
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
        let idx = thread::index_1d();
        let idx_usize = idx.get();
        if idx_usize < n as usize {
            let p = unsafe { *parents.as_mut_ptr().add(idx_usize) };
            let gp = unsafe { *parents.as_mut_ptr().add(p as usize) };
            if p != gp {
                unsafe {
                    *parents.as_mut_ptr().add(idx_usize) = gp;
                }
                let ptr = unsafe { &*(changed.as_mut_ptr() as *const DeviceAtomicU32) };
                ptr.fetch_or(1, AtomicOrdering::Relaxed);
            }
        }
    }
}
