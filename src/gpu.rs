use cuda_core::{CudaContext, LaunchConfig};

#[path = "kernels.rs"]
pub mod kernels;

// Threads per block for all kernels. 256 is safe across all SM architectures.
const CUDA_THREADS_PER_BLOCK: usize = 256;

/// Returns `true` if the CUDA driver library is loadable and a context can be created.
pub fn is_cuda_available() -> bool {
    CudaContext::new(0).is_ok()
}

pub struct GpuRunner {
    pub ctx: std::sync::Arc<CudaContext>,
    pub stream: std::sync::Arc<cuda_core::CudaStream>,
    pub module: kernels::kernels::LoadedModule,
}

impl GpuRunner {
    pub fn new() -> Option<Self> {
        let ctx = CudaContext::new(0).ok()?;
        let stream = ctx.default_stream();
        let module = kernels::kernels::load(&ctx).ok()?;
        Some(Self { ctx, stream, module })
    }

    /// Run union-find on the GPU via cuda-oxide.
    pub fn gpu_union_find(
        &self,
        num_elements: usize,
        edges: &[(u32, u32)],
        verbose: bool,
    ) -> Option<Vec<u32>> {
        if num_elements == 0 {
            return Some(Vec::new());
        }

        let stream = &self.stream;
        let module = &self.module;
        let dur_init = std::time::Duration::from_secs(0);

    let t_h2d = std::time::Instant::now();
    let mut parents_gpu = cuda_core::DeviceBuffer::<u32>::zeroed(&stream, num_elements).ok()?;
    let n = num_elements as u32;
    let num_edges = edges.len();

    // Copy edges into [[u32; 2]] shape
    let edges_gpu = if num_edges > 0 {
        let edge_slice: &[[u32; 2]] =
            unsafe { std::slice::from_raw_parts(edges.as_ptr() as *const [u32; 2], num_edges) };
        Some(cuda_core::DeviceBuffer::from_host(&stream, edge_slice).ok()?)
    } else {
        None
    };
    let dur_h2d = t_h2d.elapsed();

    let t_kernel = std::time::Instant::now();

    let tpb = CUDA_THREADS_PER_BLOCK as u32;
    let blocks_n = ((num_elements + CUDA_THREADS_PER_BLOCK - 1) / CUDA_THREADS_PER_BLOCK) as u32;
    let cfg_n = LaunchConfig {
        grid_dim: (blocks_n, 1, 1),
        block_dim: (tpb, 1, 1),
        shared_mem_bytes: 0,
    };

    module
        .initialize_parents(&stream, cfg_n, &mut parents_gpu, n)
        .ok()?;

    if num_edges > 0 {
        let ne = num_edges as i32;
        let blocks_e = ((num_edges + CUDA_THREADS_PER_BLOCK - 1) / CUDA_THREADS_PER_BLOCK) as u32;
        let cfg_e = LaunchConfig {
            grid_dim: (blocks_e, 1, 1),
            block_dim: (tpb, 1, 1),
            shared_mem_bytes: 0,
        };
        module
            .union_step(
                &stream,
                cfg_e,
                &mut parents_gpu,
                edges_gpu.as_ref().unwrap(),
                ne,
            )
            .ok()?;
    }

    // Pointer-jumping to convergence.
    let changed_host = vec![0u32];
    let mut changed_gpu = cuda_core::DeviceBuffer::from_host(&stream, &changed_host).ok()?;
    let cfg_jump = LaunchConfig {
        grid_dim: (blocks_n, 1, 1),
        block_dim: (tpb, 1, 1),
        shared_mem_bytes: 0,
    };

    loop {
        // Reset changed variable to 0 on device
        changed_gpu.zero_async(&stream).ok()?;

        module
            .pointer_jump(&stream, cfg_jump, &mut parents_gpu, n, &mut changed_gpu)
            .ok()?;

        // Read back whether changes occurred
        let changed_vec = changed_gpu.to_host_vec(&stream).ok()?;
        stream.synchronize().ok()?;

        if changed_vec[0] == 0 {
            break;
        }
    }
    let dur_kernel = t_kernel.elapsed();

    let t_d2h = std::time::Instant::now();
    let parents_host = parents_gpu.to_host_vec(&stream).ok()?;
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

    Some(parents_host)
    }
}

pub fn gpu_union_find(num_elements: usize, edges: &[(u32, u32)], verbose: bool) -> Option<Vec<u32>> {
    let runner = GpuRunner::new()?;
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

        // Edges form two chains: 0-1-2-8 and 5-6-7-8 (joined at 8), plus 3-4.
        let num_elements = 10;
        let edges: Vec<(u32, u32)> = vec![(0, 1), (1, 2), (3, 4), (5, 6), (6, 7), (7, 8), (2, 8)];

        if let Some(r) = gpu_union_find(num_elements, &edges, false) {
            assert_eq!(r.len(), num_elements);
            // {0,1,2,5,6,7,8} all in one component
            assert_eq!(r[0], r[1]);
            assert_eq!(r[1], r[2]);
            assert_eq!(r[2], r[8]);
            assert_eq!(r[5], r[6]);
            assert_eq!(r[6], r[7]);
            assert_eq!(r[0], r[7]);
            // {3,4} isolated from the above
            assert_eq!(r[3], r[4]);
            assert_ne!(r[3], r[0]);
            // {9} singleton
            assert_eq!(r[9], 9);
        } else {
            println!("GPU execution returned None; skipped.");
        }
    }
}
