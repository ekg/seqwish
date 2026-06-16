use rayon::prelude::*;
use seqwish::dset64_asm::DisjointSetsAsm;
use seqwish::gpu::gpu_union_find;
use std::time::Instant;

fn generate_mock_edges(num_elements: usize, num_edges: usize) -> Vec<(u32, u32)> {
    // Deterministic pseudo-random edge generator
    let mut state = 42u64;
    let mut next_random = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state
    };

    let mut edges = Vec::with_capacity(num_edges);
    for _ in 0..num_edges {
        let u = (next_random() % (num_elements as u64)) as u32;
        let v = (next_random() % (num_elements as u64)) as u32;
        edges.push((u, v));
    }
    edges
}

fn bench_size(num_elements: usize, num_edges: usize) {
    println!(
        "=== Benchmarking with {} elements and {} edges ===",
        num_elements, num_edges
    );
    let edges = generate_mock_edges(num_elements, num_edges);

    // CPU bench
    let start_cpu = Instant::now();
    let cpu_dsets = DisjointSetsAsm::new(num_elements);
    edges.par_iter().for_each(|&(u, v)| {
        cpu_dsets.unite(u as usize, v as usize);
    });

    // Flatten / find all roots to complete CPU disjoint-set work
    let cpu_roots: Vec<usize> = (0..num_elements)
        .into_par_iter()
        .map(|i| cpu_dsets.find(i))
        .collect();
    let duration_cpu = start_cpu.elapsed();
    println!("CPU (DisjointSetsAsm): {:?}", duration_cpu);

    // GPU bench
    let start_gpu = Instant::now();
    if let Some(gpu_roots) = gpu_union_find(num_elements, &edges, false) {
        let duration_gpu = start_gpu.elapsed();
        println!("GPU (gpu_union_find): {:?}", duration_gpu);
        println!(
            "GPU Speedup: {:.2}x",
            duration_cpu.as_secs_f64() / duration_gpu.as_secs_f64()
        );

        // Correctness verification
        let mut mismatch_count = 0;
        for i in 0..num_elements {
            let cpu_root = cpu_roots[i];
            let gpu_root = gpu_roots[i] as usize;
            if cpu_dsets.find(cpu_root) != cpu_dsets.find(gpu_root) {
                mismatch_count += 1;
            }
        }
        if mismatch_count > 0 {
            println!(
                "Warning: {} mismatching representatives found between CPU and GPU!",
                mismatch_count
            );
        } else {
            println!("Verification: CPU and GPU sets match perfectly!");
        }
    } else {
        println!("GPU: CUDA driver or device not available. Skipped GPU execution.");
    }
    println!();
}

fn main() {
    // If BENCH_CPU_THREADS is set, use that to configure Rayon.
    // Otherwise, Rayon will default to using all available logical cores on the HPC node.
    if let Ok(threads_str) = std::env::var("BENCH_CPU_THREADS") {
        if let Ok(threads) = threads_str.parse::<usize>() {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap_or(());
        }
    }

    println!("=== CPU Benchmarking Configuration ===");
    println!(
        "Rayon CPU threads in pool: {}",
        rayon::current_num_threads()
    );
    println!();

    // Small dataset
    bench_size(100_000, 200_000);

    // Medium dataset
    bench_size(1_000_000, 2_000_000);

    // Large dataset
    bench_size(10_000_000, 20_000_000);
}
