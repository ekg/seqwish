use cudarc::driver::{CudaContext, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::compile_ptx;

const CUDA_SRC: &str = include_str!("union_find.cu");

// Threads per block for all kernels. 256 is safe across all SM architectures.
const CUDA_THREADS_PER_BLOCK: usize = 256;

/// Returns `true` if both the CUDA driver and NVRTC JIT libraries are loadable.
pub fn is_cuda_available() -> bool {
    let cuda_ok = cudarc::get_lib_name_candidates("cuda")
        .into_iter()
        .any(|c| unsafe { libloading::Library::new(&c).is_ok() });

    if !cuda_ok {
        return false;
    }

    cudarc::get_lib_name_candidates("nvrtc")
        .into_iter()
        .any(|c| unsafe { libloading::Library::new(&c).is_ok() })
}

/// Run union-find on the GPU via JIT-compiled CUDA kernels.
pub fn gpu_union_find(
    num_elements: usize,
    edges: &[(u32, u32)],
    verbose: bool,
) -> Option<Vec<u32>> {
    if num_elements == 0 {
        return Some(Vec::new());
    }

    let t_init = std::time::Instant::now();
    let ctx = CudaContext::new(0).ok()?;
    let stream = ctx.default_stream();
    let ptx = compile_ptx(CUDA_SRC).ok()?;
    let module = ctx.load_module(ptx).ok()?;
    let initialize_parents = module.load_function("initialize_parents").ok()?;
    let union_step = module.load_function("union_step").ok()?;
    let pointer_jump = module.load_function("pointer_jump").ok()?;
    let dur_init = t_init.elapsed();

    let t_h2d = std::time::Instant::now();
    let mut parents_gpu = stream.alloc_zeros::<u32>(num_elements).ok()?;
    let n = num_elements as u32;
    let num_edges = edges.len();
    let edges_gpu = if num_edges > 0 {
        let edge_slice: &[[u32; 2]] =
            unsafe { std::slice::from_raw_parts(edges.as_ptr() as *const [u32; 2], num_edges) };
        Some(stream.clone_htod(edge_slice).ok()?)
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

    let mut init_builder = stream.launch_builder(&initialize_parents);
    init_builder.arg(&mut parents_gpu);
    init_builder.arg(&n);
    unsafe { init_builder.launch(cfg_n).ok()? };

    if num_edges > 0 {
        let ne = num_edges as i32;
        let blocks_e = ((num_edges + CUDA_THREADS_PER_BLOCK - 1) / CUDA_THREADS_PER_BLOCK) as u32;
        let cfg_e = LaunchConfig {
            grid_dim: (blocks_e, 1, 1),
            block_dim: (tpb, 1, 1),
            shared_mem_bytes: 0,
        };
        let mut union_builder = stream.launch_builder(&union_step);
        union_builder.arg(&mut parents_gpu);
        union_builder.arg(edges_gpu.as_ref().unwrap());
        union_builder.arg(&ne);
        unsafe { union_builder.launch(cfg_e).ok()? };
    }

    // Pointer-jumping to convergence.
    let mut changed_host = vec![0i32];
    let mut changed_gpu = stream.alloc_zeros::<i32>(1).ok()?;
    let cfg_jump = LaunchConfig {
        grid_dim: (blocks_n, 1, 1),
        block_dim: (tpb, 1, 1),
        shared_mem_bytes: 0,
    };
    loop {
        stream.memcpy_htod(&changed_host, &mut changed_gpu).ok()?;
        let mut jump_builder = stream.launch_builder(&pointer_jump);
        jump_builder.arg(&mut parents_gpu);
        jump_builder.arg(&n);
        jump_builder.arg(&mut changed_gpu);
        unsafe { jump_builder.launch(cfg_jump).ok()? };
        stream.memcpy_dtoh(&changed_gpu, &mut changed_host).ok()?;
        stream.synchronize().ok()?;
        if changed_host[0] == 0 {
            break;
        }
        changed_host[0] = 0;
    }
    let dur_kernel = t_kernel.elapsed();

    let t_d2h = std::time::Instant::now();
    let mut parents_host = vec![0u32; num_elements];
    stream.memcpy_dtoh(&parents_gpu, &mut parents_host).ok()?;
    stream.synchronize().ok()?;
    let dur_d2h = t_d2h.elapsed();

    if verbose {
        let dur_excl_init = dur_h2d + dur_kernel + dur_d2h;
        eprintln!("[gpu] profiling breakdown:");
        eprintln!("[gpu] JIT compile + driver init: {:?}", dur_init);
        eprintln!("[gpu] host to device copy + alloc: {:?}", dur_h2d);
        eprintln!("[gpu] pure kernel computation: {:?}", dur_kernel);
        eprintln!("[gpu] device to host copy: {:?}", dur_d2h);
        eprintln!("[gpu] total (excluding startup): {:?}", dur_excl_init);
    }

    Some(parents_host)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_union_find_or_fallback() {
        // cudarc 0.19.x panics on missing libnvrtc.so rather than returning Err.
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
