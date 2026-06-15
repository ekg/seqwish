// Integration test for seqwish in-memory vs disk mode
use flate2::write::GzEncoder;
use flate2::Compression;
use seqwish::intervaltree::IntervalTree;
use seqwish::*;
use std::fs::File;
use std::io::Write;
use std::sync::{Arc, Mutex, RwLock};
use std::time::{SystemTime, UNIX_EPOCH};

/// Helper function to run full pipeline and return GFA output
fn run_pipeline(in_memory: bool) -> Result<String, Box<dyn std::error::Error>> {
    // Create unique temporary paths
    let unique_id = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let fasta_path = format!("/tmp/test_pipeline_{}_{}.fasta", unique_id, in_memory);
    let paf_path = format!("/tmp/test_pipeline_{}_{}.paf.gz", unique_id, in_memory);
    let base_dir = format!("/tmp/test_pipeline_{}_{}", unique_id, in_memory);
    std::fs::create_dir_all(&base_dir)?;

    // Create test FASTA with 3 sequences
    let mut fasta_file = File::create(&fasta_path)?;
    writeln!(fasta_file, ">seq1")?;
    writeln!(fasta_file, "ACGTACGTACGT")?; // 12 bases
    writeln!(fasta_file, ">seq2")?;
    writeln!(fasta_file, "ACGTACGTACGT")?; // 12 bases (identical to seq1)
    writeln!(fasta_file, ">seq3")?;
    writeln!(fasta_file, "ACGTACGTTTTT")?; // 12 bases (8bp match, 4bp different)
    fasta_file.flush()?;

    // Build SeqIndex
    let mut seqidx = seqindex::SeqIndex::new();
    seqidx.build_index(&fasta_path)?;
    let seqidx = Arc::new(seqidx);

    // Create test PAF with alignments
    let paf_file = File::create(&paf_path)?;
    let mut gz = GzEncoder::new(paf_file, Compression::default());
    // seq1 vs seq2: perfect match
    writeln!(
        gz,
        "seq1\t12\t0\t12\t+\tseq2\t12\t0\t12\t12\t12\t60\tcg:Z:12M"
    )?;
    // seq1 vs seq3: 8bp match
    writeln!(gz, "seq1\t12\t0\t8\t+\tseq3\t12\t0\t8\t8\t8\t60\tcg:Z:8M")?;
    gz.finish()?;

    // Create alignment tree
    let mut aln_iitree_obj = if in_memory {
        intervaltree::AdaptiveTree::new_memory()?
    } else {
        let aln_path = format!("{}/aln.iitree", base_dir);
        intervaltree::AdaptiveTree::new_disk(&aln_path)?
    };
    aln_iitree_obj.open_writer()?;
    let aln_iitree = Arc::new(Mutex::new(aln_iitree_obj));

    // Process alignments
    alignments::unpack_paf_alignments(
        &paf_path,
        Arc::clone(&aln_iitree),
        Arc::clone(&seqidx),
        1,   // min_match_len
        0.0, // sparsification_factor (keep all)
        1,   // num_threads
    )?;

    // Index alignment tree
    {
        let mut tree_guard = aln_iitree.lock().unwrap();
        tree_guard.index()?;
    }

    // Create node and path trees
    let mut node_iitree_obj = if in_memory {
        intervaltree::AdaptiveTree::new_memory()?
    } else {
        let node_path = format!("{}/node.iitree", base_dir);
        intervaltree::AdaptiveTree::new_disk(&node_path)?
    };
    node_iitree_obj.open_writer()?;
    let node_iitree = Arc::new(RwLock::new(node_iitree_obj));

    let mut path_iitree_obj = if in_memory {
        intervaltree::AdaptiveTree::new_memory()?
    } else {
        let path_path = format!("{}/path.iitree", base_dir);
        intervaltree::AdaptiveTree::new_disk(&path_path)?
    };
    path_iitree_obj.open_writer()?;
    let path_iitree = Arc::new(RwLock::new(path_iitree_obj));

    // Convert Mutex to Arc for transclosure
    let aln_iitree_arc = Arc::new(
        Arc::try_unwrap(aln_iitree)
            .map_err(|_| "Cannot unwrap aln_iitree Arc")?
            .into_inner()
            .unwrap(),
    );

    // Run transclosure
    let seq_v_file = format!("{}/seq_v.txt", base_dir);
    let graph_length = transclosure::compute_transitive_closures(
        Arc::clone(&seqidx),
        aln_iitree_arc,
        &seq_v_file,
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        0,     // repeat_max
        0,     // min_repeat_dist
        10000, // transclose_batch_size
        false, // show_progress
        1,     // num_threads
    )?;

    // Index trees
    {
        let mut node_guard = node_iitree.write().unwrap();
        node_guard.index()?;
    }
    {
        let mut path_guard = path_iitree.write().unwrap();
        path_guard.index()?;
    }

    // Run compaction
    use bitvec::prelude::*;
    let mut seq_id_bv = BitVec::<u64, Lsb0>::repeat(false, graph_length + 1);
    compact::compact_nodes(
        Arc::clone(&seqidx),
        graph_length,
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &mut seq_id_bv,
        1, // num_threads
    )?;

    // Build rank/select bitvector
    let seq_id_cbv =
        links::RankSelectBitVector::from_bitvec(&seq_id_bv.iter().by_vals().collect::<Vec<bool>>());

    // Derive links
    let link_v = links::derive_links(
        Arc::clone(&seqidx),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &seq_id_cbv,
        1, // num_threads
    )?;

    // Generate GFA output
    let mut gfa_output = Vec::new();
    gfa::emit_gfa(
        &mut gfa_output,
        graph_length,
        &seq_v_file,
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &seq_id_cbv,
        Arc::clone(&seqidx),
        link_v.links(),
        1, // num_threads
    )?;

    // Cleanup
    std::fs::remove_file(&fasta_path).ok();
    std::fs::remove_file(&paf_path).ok();
    std::fs::remove_dir_all(&base_dir).ok();

    Ok(String::from_utf8(gfa_output)?)
}

#[test]
fn test_pipeline_disk_vs_memory_identical() -> Result<(), Box<dyn std::error::Error>> {
    // Run pipeline with disk-backed trees
    let disk_gfa = run_pipeline(false)?;

    // Run pipeline with in-memory trees
    let memory_gfa = run_pipeline(true)?;

    // Results should be identical
    assert_eq!(
        disk_gfa, memory_gfa,
        "Disk and in-memory modes produced different GFA outputs"
    );

    // Verify we got valid GFA output
    assert!(
        disk_gfa.starts_with("H\tVN:Z:1.0"),
        "GFA should start with header"
    );
    assert!(disk_gfa.contains("S\t"), "GFA should contain segments");
    assert!(disk_gfa.contains("P\t"), "GFA should contain paths");

    Ok(())
}

#[test]
fn test_pipeline_memory_mode_basic() -> Result<(), Box<dyn std::error::Error>> {
    // Run in-memory pipeline and verify it produces valid output
    let gfa = run_pipeline(true)?;

    // Parse GFA and verify structure
    let lines: Vec<&str> = gfa.lines().collect();

    // Count record types
    let mut num_header = 0;
    let mut num_segments = 0;
    let mut num_paths = 0;

    for line in &lines {
        if line.starts_with("H\t") {
            num_header += 1;
        } else if line.starts_with("S\t") {
            num_segments += 1;
        } else if line.starts_with("P\t") {
            num_paths += 1;
        }
    }

    // Verify we have expected structure
    assert_eq!(num_header, 1, "Should have exactly 1 header line");
    assert!(num_segments > 0, "Should have at least 1 segment");
    assert!(
        num_paths == 3,
        "Should have exactly 3 paths (one per input sequence)"
    );

    // Verify paths contain our sequence names
    let path_lines: Vec<&str> = lines
        .iter()
        .filter(|l| l.starts_with("P\t"))
        .copied()
        .collect();
    let path_str = path_lines.join("\n");
    assert!(path_str.contains("seq1"), "Should have path for seq1");
    assert!(path_str.contains("seq2"), "Should have path for seq2");
    assert!(path_str.contains("seq3"), "Should have path for seq3");

    Ok(())
}

#[test]
fn test_pipeline_disk_mode_basic() -> Result<(), Box<dyn std::error::Error>> {
    // Run disk-backed pipeline and verify it produces valid output
    let gfa = run_pipeline(false)?;

    // Parse GFA and verify structure (same checks as memory mode)
    let lines: Vec<&str> = gfa.lines().collect();

    let mut num_header = 0;
    let mut num_segments = 0;
    let mut num_paths = 0;

    for line in &lines {
        if line.starts_with("H\t") {
            num_header += 1;
        } else if line.starts_with("S\t") {
            num_segments += 1;
        } else if line.starts_with("P\t") {
            num_paths += 1;
        }
    }

    assert_eq!(num_header, 1, "Should have exactly 1 header line");
    assert!(num_segments > 0, "Should have at least 1 segment");
    assert!(num_paths == 3, "Should have exactly 3 paths");

    Ok(())
}
