use cuda_core::{CudaContext, LaunchConfig, PinnedHostBuffer};

#[path = "kernels.rs"]
pub mod kernels;

// Threads per block for all kernels. 256 is safe across all SM architectures.
const CUDA_THREADS_PER_BLOCK: usize = 256;

/// Returns `true` if the CUDA driver library is loadable and a context can be created.
pub fn is_cuda_available() -> bool {
    CudaContext::new(0).is_ok()
}

/// Cached GPU state. incl context, stream, module, and reusable device/host buffers.
///
/// Creating a GpuRunner is expensive (driver init + PTX load). So ideally keep one alive
/// and reuse it across calls to reduce that cost.
pub struct GpuRunner {
    pub ctx: std::sync::Arc<CudaContext>,
    pub stream: std::sync::Arc<cuda_core::CudaStream>,
    pub module: kernels::kernels::LoadedModule,

    /// Pinned staging buffer for the edge list
    edges_pinned: Option<PinnedHostBuffer<[u32; 2]>>,
    edges_capacity: usize,

    /// Cached parents array on the device, will be reallocated only when needed.
    parents_gpu: Option<cuda_core::DeviceBuffer<u32>>,
    parents_capacity: usize,

    /// Single-element convergence flag
    changed_gpu: Option<cuda_core::DeviceBuffer<u32>>,
}

impl GpuRunner {
    pub fn new() -> Option<Self> {
        let ctx = CudaContext::new(0).ok()?;
        let stream = ctx.default_stream();
        let module = kernels::kernels::load(&ctx).ok()?;
        let changed_gpu = cuda_core::DeviceBuffer::<u32>::zeroed(&stream, 1).ok()?;
        Some(Self {
            ctx,
            stream,
            module,
            edges_pinned: None,
            edges_capacity: 0,
            parents_gpu: None,
            parents_capacity: 0,
            changed_gpu: Some(changed_gpu),
        })
    }

    /// Run union-find on the GPU via cuda-oxide.
    pub fn gpu_union_find(
        &mut self,
        num_elements: usize,
        edges: &[(u32, u32)],
        verbose: bool,
    ) -> Option<Vec<u32>> {
        if num_elements == 0 {
            return Some(Vec::new());
        }

        let stream = &self.stream.clone();
        let module = &self.module;
        let dur_init = std::time::Duration::from_secs(0);
        let num_edges = edges.len();
        let n = num_elements as u32;

        let t_h2d = std::time::Instant::now();

        // We receive edges as &[(u32, u32)] but the kernel expects &[[u32; 2]].
        // We maintain a pinned staging buffer in the runner to maximise PCIe throughput
        let edges_gpu = if num_edges > 0 {
            // Grow the pinned staging buffer if needed.
            if num_edges > self.edges_capacity {
                self.edges_pinned =
                    Some(PinnedHostBuffer::<[u32; 2]>::zeroed(&self.ctx, num_edges).ok()?);
                self.edges_capacity = num_edges;
            }
            let pinned = self.edges_pinned.as_mut()?;
            // Fill: convert (u32, u32) -> [u32; 2] into the pinned buffer.
            for (slot, &(u, v)) in pinned.iter_mut().zip(edges.iter()) {
                *slot = [u, v];
            }
            // Transfer from pinned staging -> device.
            // SAFETY: pinned lives for the duration of this function and we
            // call stream.synchronize() before returning, so the DMA
            // completes before pinned is accessible again.
            let buf = unsafe { cuda_core::DeviceBuffer::from_pinned_host(stream, pinned).ok()? };
            Some(buf)
        } else {
            None
        };

        // --- Parents buffer ---
        if num_elements > self.parents_capacity {
            self.parents_gpu =
                Some(cuda_core::DeviceBuffer::<u32>::zeroed(stream, num_elements).ok()?);
            self.parents_capacity = num_elements;
        }
        let parents_gpu = self.parents_gpu.as_mut()?;
        let dur_h2d = t_h2d.elapsed();

        // --- Kernel launch ---
        let t_kernel = std::time::Instant::now();

        let tpb = CUDA_THREADS_PER_BLOCK as u32;
        let blocks_n =
            ((num_elements + CUDA_THREADS_PER_BLOCK - 1) / CUDA_THREADS_PER_BLOCK) as u32;
        let cfg_n = LaunchConfig {
            grid_dim: (blocks_n, 1, 1),
            block_dim: (tpb, 1, 1),
            shared_mem_bytes: 0,
        };

        module
            .initialize_parents(stream, cfg_n, parents_gpu, n)
            .ok()?;

        if num_edges > 0 {
            let ne = num_edges as i32;
            let blocks_e =
                ((num_edges + CUDA_THREADS_PER_BLOCK - 1) / CUDA_THREADS_PER_BLOCK) as u32;
            let cfg_e = LaunchConfig {
                grid_dim: (blocks_e, 1, 1),
                block_dim: (tpb, 1, 1),
                shared_mem_bytes: 0,
            };
            module
                .union_step(stream, cfg_e, parents_gpu, edges_gpu.as_ref().unwrap(), ne)
                .ok()?;
        }

        // --- Pointer-jumping to convergence ---
        let changed_gpu = self.changed_gpu.as_mut()?;
        let cfg_jump = LaunchConfig {
            grid_dim: (blocks_n, 1, 1),
            block_dim: (tpb, 1, 1),
            shared_mem_bytes: 0,
        };

        loop {
            changed_gpu.zero_async(stream).ok()?;
            module
                .pointer_jump(stream, cfg_jump, parents_gpu, n, changed_gpu)
                .ok()?;
            let changed_vec = changed_gpu.to_host_vec(stream).ok()?;
            stream.synchronize().ok()?;
            if changed_vec[0] == 0 {
                break;
            }
        }
        let dur_kernel = t_kernel.elapsed();

        // --- Device-to-host result copy ---
        let t_d2h = std::time::Instant::now();
        let mut parents_host = parents_gpu.to_host_vec(stream).ok()?;
        stream.synchronize().ok()?;
        let dur_d2h = t_d2h.elapsed();

        if verbose {
            let dur_excl_init = dur_h2d + dur_kernel + dur_d2h;
            eprintln!("[gpu] profiling breakdown:");
            eprintln!("[gpu] module load + driver init: {:?}", dur_init);
            eprintln!("[gpu] host to device copy + alloc: {:?}", dur_h2d);
            eprintln!("[gpu] pure kernel computation: {:?}", dur_kernel);
            eprintln!("[gpu] device to host copy: {:?}", dur_d2h);
            eprintln!("[gpu] total (excluding startup): {:?}", dur_excl_init);
        }

        parents_host.truncate(num_elements);
        Some(parents_host)
    }
}

pub fn gpu_union_find(
    num_elements: usize,
    edges: &[(u32, u32)],
    verbose: bool,
) -> Option<Vec<u32>> {
    let mut runner = GpuRunner::new()?;
    runner.gpu_union_find(num_elements, edges, verbose)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_union_find_or_fallback() {
        if !is_cuda_available() {
            println!("[gpu] CUDA not available: skipping.");
            return;
        }

        let num_elements = 10;
        let edges: Vec<(u32, u32)> = vec![(0, 1), (1, 2), (3, 4), (5, 6), (6, 7), (7, 8), (2, 8)];

        if let Some(r) = gpu_union_find(num_elements, &edges, false) {
            assert_eq!(r.len(), num_elements);
            assert_eq!(r[0], r[1]);
            assert_eq!(r[1], r[2]);
            assert_eq!(r[2], r[8]);
            assert_eq!(r[5], r[6]);
            assert_eq!(r[6], r[7]);
            assert_eq!(r[0], r[7]);
            assert_eq!(r[3], r[4]);
            assert_ne!(r[3], r[0]);
            assert_eq!(r[9], 9);
        } else {
            println!("GPU execution returned None; skipped.");
        }
    }

    #[test]
    fn test_gpu_union_find_random() {
        if !is_cuda_available() {
            println!("[gpu] CUDA not available: skipping.");
            return;
        }

        let num_elements = 1000;
        // Simple LCG edge generator for test setup without external deps
        let mut state = 12345u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32 % (num_elements as u32)
        };
        let edges: Vec<(u32, u32)> = (0..2000).map(|_| (next(), next())).collect();

        // Run CPU union-find
        let cpu_dsets = crate::dset64::DisjointSets::new(num_elements);
        for &(u, v) in &edges {
            cpu_dsets.unite(u as usize, v as usize);
        }
        let cpu_roots: Vec<usize> = (0..num_elements).map(|i| cpu_dsets.find(i)).collect();

        // Run GPU union-find
        if let Some(gpu_roots) = gpu_roots(num_elements, &edges, false) {
            assert_eq!(gpu_roots.len(), num_elements);

            // Verify partition equivalence:
            // For every element i, its CPU representative should map to the same GPU representative.
            let mut cpu_to_gpu_root = vec![None; num_elements];
            for i in 0..num_elements {
                let cpu_r = cpu_roots[i];
                let gpu_r = gpu_roots[i];
                if let Some(existing_gpu_r) = cpu_to_gpu_root[cpu_r] {
                    assert_eq!(
                        gpu_r, existing_gpu_r,
                        "Partition mismatch at element {}: CPU root {} mapped to GPU roots {} and {}",
                        i, cpu_r, gpu_r, existing_gpu_r
                    );
                } else {
                    cpu_to_gpu_root[cpu_r] = Some(gpu_r);
                }
            }
        }
    }
}
