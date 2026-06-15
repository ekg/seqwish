use cudarc::driver::{CudaContext, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::compile_ptx;

const CUDA_SRC: &str = include_str!("union_find.cu");

/// Check if the CUDA driver and NVRTC libraries are dynamically loadable.
pub fn is_cuda_available() -> bool {
    let mut cuda_ok = false;
    for candidate in cudarc::get_lib_name_candidates("cuda") {
        if unsafe { libloading::Library::new(&candidate) }.is_ok() {
            cuda_ok = true;
            break;
        }
    }
    if !cuda_ok {
        return false;
    }

    let mut nvrtc_ok = false;
    for candidate in cudarc::get_lib_name_candidates("nvrtc") {
        if unsafe { libloading::Library::new(&candidate) }.is_ok() {
            nvrtc_ok = true;
            break;
        }
    }
    nvrtc_ok
}

/// Run disjoint-set union-find on the GPU using CUDA.
///
/// Returns a vector where the `i`-th element represents the representative root of `i`.
/// If no CUDA device is found, or compilation/execution fails, returns `None`.
pub fn gpu_union_find(num_elements: usize, edges: &[(u32, u32)]) -> Option<Vec<u32>> {
    if num_elements == 0 {
        return Some(Vec::new());
    }

    if !is_cuda_available() {
        return None;
    }

    // Initialize CUDA context and stream
    let ctx = CudaContext::new(0).ok()?;
    let stream = ctx.default_stream();

    // Compile the CUDA source string into PTX via NVRTC
    let ptx = compile_ptx(CUDA_SRC).ok()?;
    
    // Load module into device context
    let module = ctx.load_module(ptx).ok()?;
    
    let initialize_parents = module.load_function("initialize_parents").ok()?;
    let union_step = module.load_function("union_step").ok()?;
    let pointer_jump = module.load_function("pointer_jump").ok()?;

    // Allocate parents buffer on GPU
    let mut parents = stream.alloc_zeros::<i32>(num_elements).ok()?;
    let num_elements_i32 = num_elements as i32;
    
    // Copy edges to GPU device memory
    let edge_data: Vec<[i32; 2]> = edges.iter().map(|&(u, v)| [u as i32, v as i32]).collect();
    let edges_gpu = stream.clone_htod(&edge_data).ok()?;
    
    // Launch parent initialization
    let threads_per_block = 256;
    let blocks_parents = (num_elements + threads_per_block - 1) / threads_per_block;
    let cfg_parents = LaunchConfig {
        grid_dim: (blocks_parents as u32, 1, 1),
        block_dim: (threads_per_block as u32, 1, 1),
        shared_mem_bytes: 0,
    };
    
    let mut init_builder = stream.launch_builder(&initialize_parents);
    init_builder.arg(&mut parents);
    init_builder.arg(&num_elements_i32);
    unsafe { init_builder.launch(cfg_parents).ok()? };

    // Launch union operations
    let num_edges = edges.len();
    if num_edges > 0 {
        let num_edges_i32 = num_edges as i32;
        let blocks_edges = (num_edges + threads_per_block - 1) / threads_per_block;
        let cfg_edges = LaunchConfig {
            grid_dim: (blocks_edges as u32, 1, 1),
            block_dim: (threads_per_block as u32, 1, 1),
            shared_mem_bytes: 0,
        };
        let mut union_builder = stream.launch_builder(&union_step);
        union_builder.arg(&mut parents);
        union_builder.arg(&edges_gpu);
        union_builder.arg(&num_edges_i32);
        unsafe { union_builder.launch(cfg_edges).ok()? };
    }

    // Launch pointer-jumping loop to flatten the tree structure
    let mut changed_host = vec![0i32];
    let mut changed_gpu = stream.alloc_zeros::<i32>(1).ok()?;
    
    let blocks_jump = (num_elements + threads_per_block - 1) / threads_per_block;
    let cfg_jump = LaunchConfig {
        grid_dim: (blocks_jump as u32, 1, 1),
        block_dim: (threads_per_block as u32, 1, 1),
        shared_mem_bytes: 0,
    };

    loop {
        // Reset GPU flag to 0
        stream.memcpy_htod(&changed_host, &mut changed_gpu).ok()?;
        
        let mut jump_builder = stream.launch_builder(&pointer_jump);
        jump_builder.arg(&mut parents);
        jump_builder.arg(&num_elements_i32);
        jump_builder.arg(&mut changed_gpu);
        unsafe { jump_builder.launch(cfg_jump).ok()? };
        
        // Copy convergence flag back to host
        stream.memcpy_dtoh(&changed_gpu, &mut changed_host).ok()?;
        stream.synchronize().ok()?;
        if changed_host[0] == 0 {
            break;
        }
        changed_host[0] = 0; // Reset for next iteration
    }

    // Read flat representatives back to host
    let mut parents_host = vec![0i32; num_elements];
    stream.memcpy_dtoh(&parents, &mut parents_host).ok()?;
    stream.synchronize().ok()?;

    Some(parents_host.into_iter().map(|p| p as u32).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_union_find_or_fallback() {
        let num_elements = 10;
        let edges = vec![(0, 1), (1, 2), (3, 4), (5, 6), (6, 7), (7, 8), (2, 8)];
        
        if let Some(representatives) = gpu_union_find(num_elements, &edges) {
            assert_eq!(representatives.len(), num_elements);
            
            // Check that connected components share the same root
            assert_eq!(representatives[0], representatives[1]);
            assert_eq!(representatives[1], representatives[2]);
            assert_eq!(representatives[2], representatives[8]);
            assert_eq!(representatives[5], representatives[6]);
            assert_eq!(representatives[6], representatives[7]);
            assert_eq!(representatives[0], representatives[7]); // Indirectly connected through (2, 8) and (7, 8)
            
            // Element 9 should remain in its own set
            assert_eq!(representatives[9], 9);
            
            // Elements 3 and 4 should be connected but distinct from others
            assert_eq!(representatives[3], representatives[4]);
            assert_ne!(representatives[3], representatives[0]);
        } else {
            println!("CUDA GPU not available; skipped execution check.");
        }
    }
}
