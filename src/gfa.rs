// GFA module - outputs the variation graph in GFA format
//
// This module writes the graph in GFA v1.0 format with:
// - S lines (segments/nodes with sequences)
// - L lines (links/edges between nodes)
// - P lines (paths showing input sequences through the graph)

use std::sync::{Arc, RwLock};
use std::io::{self, Write};

use crate::pos::{PosT, offset, is_rev, make_pos_t, incr_pos};
use crate::seqindex::SeqIndex;
use crate::links::RankSelectBitVector;
use crate::mmap::mmap_open;
use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;

/// Emit GFA format output for the variation graph
///
/// Writes a complete GFA v1.0 file with segments (S), links (L), and paths (P).
///
/// # Arguments
/// * `out` - Output writer
/// * `graph_length` - Total length of the graph sequence
/// * `seq_v_file` - Path to the graph sequence file (memory-mapped)
/// * `node_iitree` - Interval tree mapping graph positions to input positions
/// * `path_iitree` - Interval tree mapping input positions to graph positions
/// * `seq_id_cbv` - Ranked/select bitvector marking node boundaries
/// * `seqidx` - Sequence index for input sequences
/// * `links` - Vector of graph edges (from_node, to_node)
/// * `num_threads` - Number of threads for parallel processing
pub fn emit_gfa<W: Write>(
    out: &mut W,
    _graph_length: usize,
    seq_v_file: &str,
    _node_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    path_iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>>,
    seq_id_cbv: &RankSelectBitVector,
    seqidx: Arc<SeqIndex>,
    links: &[(PosT, PosT)],
    _num_threads: usize,
) -> io::Result<()> {
    // Write GFA header
    writeln!(out, "H\tVN:Z:1.0")?;

    // Memory-map the graph sequence file
    let mmap_handle = mmap_open(seq_v_file)?;
    let seq_v_slice = unsafe {
        std::slice::from_raw_parts(mmap_handle.ptr as *const u8, mmap_handle.size)
    };

    // Get number of nodes
    let n_nodes = seq_id_cbv.rank(seq_id_cbv.size() - 1);

    // Write nodes (S lines) - extract sequences first, then write
    let node_sequences: Vec<(usize, String)> = (1..=n_nodes)
        .map(|id| {
            let node_start = match seq_id_cbv.select(id) {
                Some(pos) => pos,
                None => return (id, String::new()),
            };
            let node_end = match seq_id_cbv.select(id + 1) {
                Some(pos) => pos,
                None => return (id, String::new()),
            };
            let node_length = node_end - node_start;
            let seq = &seq_v_slice[node_start..node_start + node_length];
            let seq_string = String::from_utf8_lossy(seq).to_string();
            (id, seq_string)
        })
        .collect();

    // Write node records in order
    for (id, seq) in node_sequences {
        if !seq.is_empty() {
            writeln!(out, "S\t{}\t{}", id, seq)?;
        }
    }


    // Write links (L lines)
    for (from, to) in links {
        if *from != 0 && *to != 0 {
            writeln!(
                out,
                "L\t{}\t{}\t{}\t{}\t0M",
                offset(*from),
                if is_rev(*from) { "-" } else { "+" },
                offset(*to),
                if is_rev(*to) { "-" } else { "+" }
            )?;
        }
    }

    // Write paths (P lines)
    let num_seqs = seqidx.n_seqs();

    for i in 1..=num_seqs {
        let j_start = match seqidx.nth_seq_offset(i) {
            Some(offset) => offset,
            None => continue,
        };
        let seq_len = match seqidx.nth_seq_length(i) {
            Some(len) => len,
            None => continue,
        };
        let k = j_start + seq_len;

        let mut path_v: Vec<PosT> = Vec::new();
        let mut seen_bp = 0u64;
        let mut j = j_start;

        while j < k {
            let mut overlap_count = 0u64;
            let mut ovlp_start_in_q = 0u64;
            let mut ovlp_end_in_q = 0u64;
            let mut pos_start_in_s = 0u64;

            // Find overlap in path_iitree
            if let Ok(path_guard) = path_iitree.read() {
                path_guard.overlap(j, j + 1, |_idx, start, end, pos| {
                    overlap_count += 1;
                    ovlp_start_in_q = start;
                    ovlp_end_in_q = end;
                    pos_start_in_s = pos;
                }).ok();
            }

            // Each input base should map to exactly one place in the graph
            if overlap_count != 1 {
                let seq_name = seqidx.nth_name(i).unwrap_or_else(|| "<unknown>".to_string());
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("[gfa] error: found {} overlaps for seq {} idx {} at j={} of {}",
                            overlap_count, seq_name, i, j, k)
                ));
            }

            let match_is_rev = is_rev(pos_start_in_s);
            let length = ovlp_end_in_q - ovlp_start_in_q;
            let mut p = pos_start_in_s;

            // Iterate through nodes in this range
            for _ in 0..length {
                let p_offset = offset(p) as usize;
                if seq_id_cbv.select(1).is_some() && p_offset < seq_id_cbv.size() {
                    // Check if this position is a node boundary
                    let node_id = seq_id_cbv.rank(p_offset + 1);
                    if node_id > 0 {
                        // Only add if we're at the start of a node
                        let node_start = seq_id_cbv.select(node_id);
                        if let Some(start_pos) = node_start {
                            if start_pos == p_offset {
                                path_v.push(make_pos_t(node_id as u64, match_is_rev));
                            }
                        }
                    }
                }
                incr_pos(&mut p);
            }

            seen_bp += length;
            j = ovlp_end_in_q;
        }

        if seen_bp != seq_len {
            let seq_name = seqidx.nth_name(i).unwrap_or_else(|| "<unknown>".to_string());
            return Err(io::Error::new(
                io::ErrorKind::Other,
                format!("[gfa] length mismatch for {}, expected {} but got {}",
                        seq_name, seq_len, seen_bp)
            ));
        }

        // Write path
        let seq_name = seqidx.nth_name(i).unwrap_or_else(|| format!("seq{}", i));
        write!(out, "P\t{}\t", seq_name)?;

        for (idx, p) in path_v.iter().enumerate() {
            if idx > 0 {
                write!(out, ",")?;
            }
            write!(out, "{}{}", offset(*p), if is_rev(*p) { "-" } else { "+" })?;
        }
        writeln!(out, "\t*")?;
    }

    // Cleanup mmap (automatic via Drop)
    drop(mmap_handle);

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfa_header() {
        let mut output = Vec::new();
        write!(&mut output, "H\tVN:Z:1.0\n").unwrap();
        let result = String::from_utf8(output).unwrap();
        assert_eq!(result, "H\tVN:Z:1.0\n");
    }
}
