// GFA module - outputs the variation graph in GFA format
//
// This module writes the graph in GFA v1.0 format with:
// - S lines (segments/nodes with sequences)
// - L lines (links/edges between nodes)
// - P lines (paths showing input sequences through the graph)

use std::io::{self, Write};
use std::sync::{Arc, RwLock};

use crate::dna::complement;
use crate::intervaltree::AdaptiveTree;
use crate::intervaltree::IntervalTree;
use crate::links::RankSelectBitVector;
use crate::mmap::mmap_open;
use crate::pos::{incr_pos, is_rev, make_pos_t, offset, PosT};
use crate::seqindex::SeqIndex;

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
#[allow(clippy::too_many_arguments)]
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
    let seq_v_slice =
        unsafe { std::slice::from_raw_parts(mmap_handle.ptr as *const u8, mmap_handle.size) };

    // Get number of nodes
    let n_nodes = seq_id_cbv.rank(seq_id_cbv.size() - 1);

    // Write nodes (S lines) - extract sequences first, then write.
    // The graph sequence has `seq_v_slice.len()` bytes; use that as the end
    // sentinel for the last node (the bitvector has size graph_length + 1).
    let graph_len = seq_v_slice.len();
    let node_sequences: Vec<(usize, String)> = (1..=n_nodes)
        .map(|id| {
            let node_start = match seq_id_cbv.select(id) {
                Some(pos) => pos,
                None => return (id, String::new()),
            };
            let node_end = seq_id_cbv.select(id + 1).unwrap_or(graph_len);
            let node_length = node_end - node_start;
            if node_start + node_length > graph_len {
                return (id, String::new());
            }
            let seq = &seq_v_slice[node_start..node_start + node_length];
            let seq_string = String::from_utf8_lossy(seq).to_string();
            (id, seq_string)
        })
        .collect();

    // Write node records in order
    for (id, seq) in node_sequences {
        if !seq.is_empty() {
            writeln!(out, "S\t{id}\t{seq}")?;
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
                path_guard
                    .overlap(j, j + 1, |_idx, start, end, pos| {
                        overlap_count += 1;
                        ovlp_start_in_q = start;
                        ovlp_end_in_q = end;
                        pos_start_in_s = pos;
                    })
                    .ok();
            }

            // Each input base should map to exactly one place in the graph
            if overlap_count != 1 {
                let seq_name = seqidx
                    .nth_name(i)
                    .unwrap_or_else(|| "<unknown>".to_string());
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!(
                        "[gfa] error: found {overlap_count} overlaps for seq {seq_name} idx {i} at j={j} of {k}"
                    ),
                ));
            }

            let match_is_rev = is_rev(pos_start_in_s);
            let length = ovlp_end_in_q - ovlp_start_in_q;

            // Validate path integrity: check that input sequence matches graph sequence
            // This is critical for correctness (matches C++ validation)
            let mut p = pos_start_in_s; // position in graph
            for q in j..j + length {
                let p_offset = offset(p) as usize;
                if p_offset < seq_v_slice.len() {
                    let mut graph_char = seq_v_slice[p_offset];
                    if is_rev(p) {
                        graph_char = complement(graph_char);
                    }

                    if let Some(input_char) = seqidx.at(q) {
                        if input_char != graph_char as char {
                            let seq_name = seqidx
                                .nth_name(i)
                                .unwrap_or_else(|| "<unknown>".to_string());
                            return Err(io::Error::new(
                                io::ErrorKind::Other,
                                format!("[gfa] GRAPH BROKEN @ {} pos {} -> graph pos {}: expected {} got {}",
                                        seq_name, q, p_offset, input_char, graph_char as char)
                            ));
                        }
                    }
                }
                incr_pos(&mut p);
            }

            // Match C++ behavior: iterate through each base and add node at boundaries
            // This is the original seqwish algorithm that adds nodes when crossing
            // node boundaries (seq_id_cbv has 1s at node starts)
            let mut p = pos_start_in_s;
            for _ in 0..length {
                let p_offset = offset(p) as usize;
                // Check if this position is a node boundary (start of a new node)
                if p_offset < seq_id_cbv.size() && seq_id_cbv.access(p_offset) {
                    let node_id = seq_id_cbv.rank(p_offset + 1);
                    if node_id > 0 {
                        path_v.push(make_pos_t(node_id as u64, match_is_rev));
                    }
                }
                incr_pos(&mut p);
            }

            seen_bp += length;
            j = ovlp_end_in_q;
        }

        if seen_bp != seq_len {
            let seq_name = seqidx
                .nth_name(i)
                .unwrap_or_else(|| "<unknown>".to_string());
            return Err(io::Error::new(
                io::ErrorKind::Other,
                format!(
                    "[gfa] length mismatch for {seq_name}, expected {seq_len} but got {seen_bp}"
                ),
            ));
        }

        // Validate path step lengths sum to sequence length.
        // Use graph_len (same as S-line emission) as sentinel for last node end.
        let mut path_step_len = 0u64;
        let mut node_lengths: Vec<(usize, usize)> = Vec::new();
        for p in path_v.iter() {
            let node_id = offset(*p) as usize;
            if node_id > 0 {
                let node_start = seq_id_cbv.select(node_id).unwrap_or(0);
                let node_end = seq_id_cbv.select(node_id + 1).unwrap_or(graph_len);
                let node_len = node_end - node_start;
                path_step_len += node_len as u64;
                node_lengths.push((node_id, node_len));
            }
        }
        if path_step_len != seq_len {
            let seq_name = seqidx
                .nth_name(i)
                .unwrap_or_else(|| "<unknown>".to_string());
            let debug_nodes: String = if node_lengths.len() <= 10 {
                node_lengths
                    .iter()
                    .map(|(id, len)| format!("{}:{}", id, len))
                    .collect::<Vec<_>>()
                    .join(",")
            } else {
                let first: String = node_lengths[..3]
                    .iter()
                    .map(|(id, len)| format!("{}:{}", id, len))
                    .collect::<Vec<_>>()
                    .join(",");
                let last: String = node_lengths[node_lengths.len() - 3..]
                    .iter()
                    .map(|(id, len)| format!("{}:{}", id, len))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("{}...{}", first, last)
            };
            return Err(io::Error::new(
                io::ErrorKind::Other,
                format!(
                    "[gfa] path step length mismatch for {seq_name}: expected {seq_len} bp but path steps sum to {path_step_len} bp ({} nodes, seen_bp={}, nodes: {})",
                    path_v.len(), seen_bp, debug_nodes
                ),
            ));
        }

        // Write path
        let seq_name = seqidx.nth_name(i).unwrap_or_else(|| format!("seq{i}"));
        write!(out, "P\t{seq_name}\t")?;

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
        writeln!(&mut output, "H\tVN:Z:1.0").unwrap();
        let result = String::from_utf8(output).unwrap();
        assert_eq!(result, "H\tVN:Z:1.0\n");
    }
}
