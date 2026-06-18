use rayon::prelude::*;
use seqwish::dset64_asm::DisjointSetsAsm;
use seqwish::gpu::GpuRunner;
use std::time::Instant;

// ---------------------------------------------------------------------------
// Non-atomic sequential union-find (fair single-threaded CPU baseline)
// Uses plain array reads/writes with no atomics.
// ---------------------------------------------------------------------------
struct PlainUnionFind {
    parent: Vec<u32>,
}

impl PlainUnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            // Path-halving: point to grandparent
            let gp = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = gp;
            x = gp;
        }
        x
    }

    fn unite(&mut self, a: u32, b: u32) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        // Union by index (lower root wins)
        if ra < rb {
            self.parent[ra as usize] = rb;
        } else {
            self.parent[rb as usize] = ra;
        }
    }

    fn roots(&mut self, n: usize) -> Vec<u32> {
        (0..n as u32).map(|i| self.find(i)).collect()
    }
}

// ---------------------------------------------------------------------------
// Edge generators
// ---------------------------------------------------------------------------

/// LCG-based generator: produces numerically correlated edges.
/// Might favour CPU cache prefetcher.
fn generate_lcg_edges(num_elements: usize, num_edges: usize) -> Vec<(u32, u32)> {
    let mut state = 42u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u32 % (num_elements as u32)
    };
    (0..num_edges).map(|_| (next(), next())).collect()
}

/// Uniform random generator using splitmix64.
/// Hopefully more representative of long-range pangenome alignments.
fn generate_uniform_edges(num_elements: usize, num_edges: usize) -> Vec<(u32, u32)> {
    let mut state = 0x9e3779b97f4a7c15u64;
    let mut next = move || {
        state = state.wrapping_add(0x9e3779b97f4a7c15);
        let mut x = state;
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
        (x ^ (x >> 31)) as u32 % (num_elements as u32)
    };
    (0..num_edges).map(|_| (next(), next())).collect()
}

// ---------------------------------------------------------------------------
// Benchmark for one edge set
// ---------------------------------------------------------------------------
fn bench_one(
    runner: Option<&mut GpuRunner>,
    label: &str,
    num_elements: usize,
    edges: &[(u32, u32)],
) {
    println!("  [{label}]");

    // --- Non-atomic sequential ---
    let t = Instant::now();
    let mut uf_seq = PlainUnionFind::new(num_elements);
    for &(u, v) in edges {
        uf_seq.unite(u, v);
    }
    let seq_roots = uf_seq.roots(num_elements);
    let dur_seq = t.elapsed();
    println!("  CPU sequential (plain, no atomics): {:?}", dur_seq);

    // --- Parallel with atomics (DisjointSetsAsm) ---
    let t = Instant::now();
    let dsets_par = DisjointSetsAsm::new(num_elements);
    edges.par_iter().for_each(|&(u, v)| {
        dsets_par.unite(u as usize, v as usize);
    });
    let _par_roots: Vec<usize> = (0..num_elements)
        .into_par_iter()
        .map(|i| dsets_par.find(i))
        .collect();
    let dur_par = t.elapsed();
    println!(
        "  CPU parallel (atomic, {} threads): {:?}",
        rayon::current_num_threads(),
        dur_par
    );

    // --- GPU ---
    let t = Instant::now();
    let gpu_roots_opt = if let Some(r) = runner {
        r.gpu_union_find(num_elements, edges, true)
    } else {
        None
    };
    if let Some(gpu_roots) = gpu_roots_opt {
        let dur_gpu = t.elapsed();
        println!("  GPU total (incl. PCIe transfers): {:?}", dur_gpu);
        println!(
            "  GPU speedup vs sequential (fair):  {:.2}x",
            dur_seq.as_secs_f64() / dur_gpu.as_secs_f64()
        );
        println!(
            "  GPU speedup vs parallel (atomic):  {:.2}x",
            dur_par.as_secs_f64() / dur_gpu.as_secs_f64()
        );

        // Correctness: check that each element lands in the same component
        let mut mismatches = 0usize;
        let mut uf_check = PlainUnionFind::new(num_elements);
        for &(u, v) in edges {
            uf_check.unite(u, v);
        }
        for i in 0..num_elements {
            let cpu_root = uf_check.find(seq_roots[i]) as usize;
            let gpu_root = uf_check.find(gpu_roots[i]) as usize;
            if cpu_root != gpu_root {
                mismatches += 1;
            }
        }
        if mismatches == 0 {
            println!("  Verification: ✓ CPU and GPU sets match");
        } else {
            println!("  Verification: ✗ {} mismatches!", mismatches);
        }
    } else {
        println!("  GPU: not available — skipped.");
    }
}

// ---------------------------------------------------------------------------
// Benchmark one dataset size with both edge distributions
// ---------------------------------------------------------------------------
fn bench_size(runner: Option<&mut GpuRunner>, num_elements: usize, num_edges: usize) {
    println!("\n=== {} elements, {} edges ===", num_elements, num_edges);

    let lcg_edges = generate_lcg_edges(num_elements, num_edges);
    let uni_edges = generate_uniform_edges(num_elements, num_edges);

    // Use a pointer shenanigan to re-borrow. Safe because the calls are sequential.
    if let Some(r) = runner {
        bench_one(
            Some(r),
            "LCG / clustered edges (CPU-friendly)",
            num_elements,
            &lcg_edges,
        );
        bench_one(
            Some(r),
            "Uniform random edges (realistic)",
            num_elements,
            &uni_edges,
        );
    } else {
        bench_one(
            None,
            "LCG / clustered edges (CPU-friendly)",
            num_elements,
            &lcg_edges,
        );
        bench_one(
            None,
            "Uniform random edges (realistic)",
            num_elements,
            &uni_edges,
        );
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
fn main() {
    // Respect BENCH_CPU_THREADS if set
    if let Ok(s) = std::env::var("BENCH_CPU_THREADS") {
        if let Ok(n) = s.parse::<usize>() {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build_global()
                .unwrap_or(());
        }
    }

    println!(
        "=== CPU: {} Rayon threads ===",
        rayon::current_num_threads()
    );

    let mut runner = GpuRunner::new();
    if runner.is_none() {
        println!("WARNING: GPU not available — GPU benchmarks will be skipped.");
    }

    // GPU warm-up: initialise driver, load PTX, and pre-allocate buffers
    // at a minimal size. The cached buffers should grow on the first real call.
    if let Some(r) = &mut runner {
        r.gpu_union_find(10, &[(0, 1)], false);
    }

    bench_size(runner.as_mut(), 100_000, 200_000);
    bench_size(runner.as_mut(), 1_000_000, 2_000_000);
    bench_size(runner.as_mut(), 10_000_000, 20_000_000);
}
