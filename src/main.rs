// seqwish: a variation graph inducer
//
// Complete Rust implementation of the seqwish algorithm for building
// variation graphs from pairwise alignments.

// Allow dead code - some functions are exported for FFI or kept for future use
#![allow(dead_code)]
// Allow complex types and many arguments in certain contexts
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]

use std::fs::File;
use std::io::{self, BufWriter};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use bitvec::prelude::*;
use clap::{Arg, Command};

mod alignments;
mod cigar;
mod compact;
mod dna;
mod dset64;
mod dset64_asm;
mod dset64_unsafe;
mod gfa;
#[cfg(feature = "cuda")]
mod gpu;

#[cfg(not(feature = "cuda"))]
mod gpu {
    pub fn is_cuda_available() -> bool {
        false
    }
    pub fn gpu_union_find(
        _num_elements: usize,
        _edges: &[(u32, u32)],
        _verbose: bool,
    ) -> Option<Vec<u32>> {
        None
    }
    pub struct GpuRunner;
    impl GpuRunner {
        pub fn new() -> Option<Self> {
            None
        }
        pub fn gpu_union_find(
            &mut self,
            _n: usize,
            _e: &[(u32, u32)],
            _v: bool,
        ) -> Option<Vec<u32>> {
            None
        }
    }
}
mod intervaltree;
mod links;
mod mmap;
mod paf;
mod pos;
mod seqindex;
mod sxs;
mod tempfile;
mod time;
mod transclosure;
mod utils;
mod version;

use crate::alignments::unpack_paf_alignments;
use crate::compact::compact_nodes;
use crate::gfa::emit_gfa;
use crate::intervaltree::{AdaptiveTree, IntervalTree};
use crate::links::{derive_links, RankSelectBitVector};
use crate::seqindex::SeqIndex;
use crate::transclosure::compute_transitive_closures;

fn main() -> io::Result<()> {
    let matches = Command::new("seqwish")
        .version(version::get_version())
        .about("A variation graph inducer")
        .arg(Arg::new("paf-alns")
            .short('p')
            .long("paf-alns")
            .value_name("FILE")
            .help("Induce the graph from these PAF formatted alignments")
            .required(true))
        .arg(Arg::new("seqs")
            .short('s')
            .long("seqs")
            .value_name("FILE")
            .help("The sequences used to generate the alignments (FASTA, FASTQ, .seq)")
            .required(true))
        .arg(Arg::new("gfa")
            .short('g')
            .long("gfa")
            .value_name("FILE")
            .help("Write the graph in GFA to FILE (stdout if not specified)"))
        .arg(Arg::new("temp-dir")
            .short('b')
            .long("temp-dir")
            .value_name("PATH")
            .help("Directory for temporary files [default: current directory]"))
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .value_name("N")
            .help("Use this many threads during parallel steps")
            .default_value("1"))
        .arg(Arg::new("repeat-max")
            .short('r')
            .long("repeat-max")
            .value_name("N")
            .help("Limit transitive closure to include no more than N copies of a given input base")
            .default_value("0"))
        .arg(Arg::new("min-repeat-distance")
            .short('l')
            .long("min-repeat-distance")
            .value_name("N")
            .help("Prevent transitive closure for bases at least this far apart in input sequences")
            .default_value("0"))
        .arg(Arg::new("min-match-len")
            .short('k')
            .long("min-match-len")
            .value_name("N")
            .help("Filter exact matches below this length")
            .default_value("0"))
        .arg(Arg::new("sparse-factor")
            .short('f')
            .long("sparse-factor")
            .value_name("N")
            .help("Sparsify input matches, keeping the fraction that minimize a hash function")
            .default_value("0.0"))
        .arg(Arg::new("transclose-batch")
            .short('B')
            .long("transclose-batch")
            .value_name("N")
            .help("Number of bp to use for transitive closure batch")
            .default_value("1000000"))
        .arg(Arg::new("keep-temp")
            .short('T')
            .long("keep-temp")
            .help("Keep intermediate files generated during graph induction")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("show-progress")
            .short('P')
            .long("show-progress")
            .help("Log algorithm progress")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("in-memory")
            .short('M')
            .long("in-memory")
            .help("Use in-memory interval trees instead of disk-backed (faster for small datasets)")
            .action(clap::ArgAction::SetTrue))
        .get_matches();

    let start_time = Instant::now();

    // Parse arguments
    let paf_file = matches.get_one::<String>("paf-alns").unwrap();
    let seq_file = matches.get_one::<String>("seqs").unwrap();
    let gfa_file = matches.get_one::<String>("gfa");
    let num_threads: usize = matches
        .get_one::<String>("threads")
        .unwrap()
        .parse()
        .unwrap_or(1);

    // Configure Rayon's global thread pool to respect num_threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    let repeat_max: u64 = matches
        .get_one::<String>("repeat-max")
        .unwrap()
        .parse()
        .unwrap_or(0);
    let min_repeat_dist: u64 = matches
        .get_one::<String>("min-repeat-distance")
        .unwrap()
        .parse()
        .unwrap_or(0);
    let min_match_len: u64 = matches
        .get_one::<String>("min-match-len")
        .unwrap()
        .parse()
        .unwrap_or(0);
    let sparse_factor: f32 = matches
        .get_one::<String>("sparse-factor")
        .unwrap()
        .parse()
        .unwrap_or(0.0);
    let transclose_batch: u64 = matches
        .get_one::<String>("transclose-batch")
        .unwrap()
        .parse()
        .unwrap_or(1000000);
    let keep_temp = matches.get_flag("keep-temp");
    let show_progress = matches.get_flag("show-progress");
    let use_in_memory = matches.get_flag("in-memory");

    // Set up temp directory
    if let Some(temp_dir) = matches.get_one::<String>("temp-dir") {
        tempfile::set_dir(temp_dir);
    }
    tempfile::set_keep_temp(keep_temp);

    // Check input files exist
    if !std::path::Path::new(seq_file).exists() {
        eprintln!("[seqwish] ERROR: input sequence file {seq_file} does not exist");
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "sequence file not found",
        ));
    }
    if !std::path::Path::new(paf_file).exists() {
        eprintln!("[seqwish] ERROR: input alignment file {paf_file} does not exist");
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "alignment file not found",
        ));
    }

    // 1) Build the sequence index
    if show_progress {
        eprintln!(
            "[seqwish::seqindex] {:.3} loading sequences",
            start_time.elapsed().as_secs_f64()
        );
    }
    let mut seqidx = SeqIndex::new();
    seqidx
        .build_index(seq_file)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    let seqidx = Arc::new(seqidx);
    if show_progress {
        eprintln!(
            "[seqwish::seqindex] {:.3} loaded {} sequences",
            start_time.elapsed().as_secs_f64(),
            seqidx.n_seqs()
        );
    }

    // 2) Index alignments
    if show_progress {
        eprintln!(
            "[seqwish::alignments] {:.3} loading alignments",
            start_time.elapsed().as_secs_f64()
        );
    }
    let aln_iitree_idx = tempfile::create("seqwish-", ".sqa")?;
    let mut aln_iitree_obj = if use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&aln_iitree_idx)?
    };
    aln_iitree_obj.open_writer()?;
    let aln_iitree = Arc::new(Mutex::new(aln_iitree_obj));

    unpack_paf_alignments(
        paf_file,
        Arc::clone(&aln_iitree),
        Arc::clone(&seqidx),
        min_match_len,
        sparse_factor,
        num_threads,
    )?;

    if show_progress {
        eprintln!(
            "[seqwish::alignments] {:.3} indexing",
            start_time.elapsed().as_secs_f64()
        );
    }
    aln_iitree.lock().unwrap().index()?;
    if show_progress {
        eprintln!(
            "[seqwish::alignments] {:.3} index built",
            start_time.elapsed().as_secs_f64()
        );
    }

    // Unwrap the Mutex - aln_iitree is now read-only, no need for mutex
    let aln_iitree_readonly = Arc::new(
        Arc::try_unwrap(aln_iitree)
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Failed to unwrap Arc"))?
            .into_inner()
            .unwrap(),
    );

    // 3) Find transitive closures and construct graph sequence
    if show_progress {
        eprintln!(
            "[seqwish::transclosure] {:.3} computing transitive closures",
            start_time.elapsed().as_secs_f64()
        );
    }
    let seq_v_file = tempfile::create("seqwish-", ".sqs")?;
    let node_iitree_idx = tempfile::create("seqwish-", ".sqn")?;
    let path_iitree_idx = tempfile::create("seqwish-", ".sqp")?;

    let mut node_iitree_obj = if use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&node_iitree_idx)?
    };
    node_iitree_obj.open_writer()?;
    let node_iitree = Arc::new(std::sync::RwLock::new(node_iitree_obj));

    let mut path_iitree_obj = if use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&path_iitree_idx)?
    };
    path_iitree_obj.open_writer()?;
    let path_iitree = Arc::new(std::sync::RwLock::new(path_iitree_obj));

    let graph_length = compute_transitive_closures(
        Arc::clone(&seqidx),
        Arc::clone(&aln_iitree_readonly),
        seq_v_file.to_str().unwrap(),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        repeat_max,
        min_repeat_dist,
        transclose_batch,
        show_progress,
        num_threads,
    )?;

    if show_progress {
        eprintln!(
            "[seqwish::transclosure] {:.3} done with transitive closures (graph length: {})",
            start_time.elapsed().as_secs_f64(),
            graph_length
        );
    }

    // 4) Compact nodes by marking boundaries
    if show_progress {
        eprintln!(
            "[seqwish::compact] {:.3} compacting nodes",
            start_time.elapsed().as_secs_f64()
        );
    }
    let mut seq_id_bv = BitVec::<u64, Lsb0>::repeat(false, graph_length + 1);

    compact_nodes(
        Arc::clone(&seqidx),
        graph_length,
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &mut seq_id_bv,
        num_threads,
    )?;

    if show_progress {
        eprintln!(
            "[seqwish::compact] {:.3} done compacting",
            start_time.elapsed().as_secs_f64()
        );
    }

    // Build rank/select structure
    let seq_id_cbv =
        RankSelectBitVector::from_bitvec(&seq_id_bv.iter().by_vals().collect::<Vec<bool>>());
    drop(seq_id_bv); // Free memory

    if show_progress {
        eprintln!(
            "[seqwish::compact] {:.3} built node index",
            start_time.elapsed().as_secs_f64()
        );
    }

    // 5) Derive links between nodes
    if show_progress {
        eprintln!(
            "[seqwish::links] {:.3} finding graph links",
            start_time.elapsed().as_secs_f64()
        );
    }

    let link_set = derive_links(
        Arc::clone(&seqidx),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &seq_id_cbv,
        num_threads,
    )?;

    if show_progress {
        eprintln!(
            "[seqwish::links] {:.3} links derived ({} links)",
            start_time.elapsed().as_secs_f64(),
            link_set.len()
        );
    }

    // 6) Emit the graph in GFA format
    if show_progress {
        eprintln!(
            "[seqwish::gfa] {:.3} writing graph",
            start_time.elapsed().as_secs_f64()
        );
    }

    if let Some(gfa_path) = gfa_file {
        // Use large buffer (1MB) for GFA output to minimize write syscalls
        let mut out = BufWriter::with_capacity(1024 * 1024, File::create(gfa_path)?);
        emit_gfa(
            &mut out,
            graph_length,
            seq_v_file.to_str().unwrap(),
            Arc::clone(&node_iitree),
            Arc::clone(&path_iitree),
            &seq_id_cbv,
            Arc::clone(&seqidx),
            link_set.links(),
            num_threads,
        )?;
    } else {
        let stdout = io::stdout();
        // Use large buffer (1MB) for stdout GFA output to minimize write syscalls
        let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
        emit_gfa(
            &mut out,
            graph_length,
            seq_v_file.to_str().unwrap(),
            Arc::clone(&node_iitree),
            Arc::clone(&path_iitree),
            &seq_id_cbv,
            Arc::clone(&seqidx),
            link_set.links(),
            num_threads,
        )?;
    }

    if show_progress {
        eprintln!(
            "[seqwish::gfa] {:.3} done",
            start_time.elapsed().as_secs_f64()
        );
    }

    Ok(())
}
