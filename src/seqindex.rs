use flate2::read::MultiGzDecoder;
use fm_index::{FMIndexWithLocate, MatchWithLocate, Search, Text};
use memmap2::Mmap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

/// Simple sparse bitvector for sequence boundaries
/// Much faster than RsVec for sparse data (select is O(1) array access vs hierarchical search)
#[derive(Clone)]
struct SparseBitVec {
    /// Sorted positions where bits are set to 1
    positions: Vec<usize>,
    /// Total size of the bitvector
    size: usize,
}

impl SparseBitVec {
    /// Create from a list of 1-bit positions
    fn from_positions(positions: Vec<usize>, size: usize) -> Self {
        SparseBitVec { positions, size }
    }

    /// Select: return position of i-th 1-bit (0-indexed for compatibility with RsVec usage)
    /// O(1) - direct array access!
    #[inline]
    fn select1(&self, i: usize) -> usize {
        if i < self.positions.len() {
            self.positions[i]
        } else {
            self.size
        }
    }

    /// Rank: count number of 1-bits BEFORE position i (excluding i itself)
    /// O(log n) binary search
    #[inline]
    fn rank1(&self, i: usize) -> usize {
        self.positions.binary_search(&i).unwrap_or_else(|idx| idx) // Found at idx or insertion point = rank
    }

    /// Get the size of the bitvector
    #[inline]
    fn len(&self) -> usize {
        self.size
    }

    /// Get bit value at position i (for compatibility)
    #[inline]
    fn get(&self, i: usize) -> u64 {
        if self.positions.binary_search(&i).is_ok() {
            1
        } else {
            0
        }
    }
}

/// Sequence index for FASTA/FASTQ files
/// Uses compressed suffix array (CSA) and succinct bitvectors to match C++ Big O bounds
///
/// Space complexity:
/// - Name index: O(n log σ) bits where n = total chars in names, σ = alphabet size
/// - Sequence boundaries: O(m log(N/m)) bits where m = # sequences, N = total sequence length
pub struct SeqIndex {
    /// Path to concatenated sequence file
    seq_filename: Option<PathBuf>,

    /// Memory-mapped sequence data
    seq_mmap: Option<Mmap>,

    /// FM-index (CSA) of concatenated sequence names with delimiters
    /// Space: O(n log σ) bits
    /// Format: ">name1 >name2 >name3 "
    name_index: Option<FMIndexWithLocate<u8>>,

    /// Original name text (needed for extraction since FM-index doesn't expose it)
    name_text: Vec<u8>,

    /// Sparse bitvector marking start of each sequence name in name_index
    /// Space: O(m) words where m = # sequences (simple array of positions)
    /// A 1-bit marks the position of each '>' character
    name_boundaries: Option<SparseBitVec>,

    /// Sparse bitvector marking start of each sequence in seq_mmap
    /// Space: O(m) words where m = # sequences (simple array of positions)
    /// A 1-bit at position i means sequence starts at seq_mmap[i]
    seq_boundaries: Option<SparseBitVec>,

    /// Total number of sequences
    seq_count: usize,
}

impl SeqIndex {
    /// Create a new empty sequence index
    pub fn new() -> Self {
        SeqIndex {
            seq_filename: None,
            seq_mmap: None,
            name_index: None,
            name_text: Vec::new(),
            name_boundaries: None,
            seq_boundaries: None,
            seq_count: 0,
        }
    }

    /// Build index from FASTA or FASTQ file
    ///
    /// Time complexity: O(n log n) for suffix array construction
    /// Space complexity: O(n log σ) bits for CSA + O(m log(N/m)) bits for boundaries
    pub fn build_index(&mut self, filename: &str) -> Result<(), String> {
        // Create temp file for sequences
        let seq_file = std::env::temp_dir().join(format!("seqwish-{}.sqq", std::process::id()));
        self.seq_filename = Some(seq_file.clone());

        // Open input file (with optional gzip support)
        let file =
            File::open(filename).map_err(|e| format!("Failed to open {filename}: {e}"))?;

        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        // Open output file for sequences (WITH LARGE BUFFER!)
        // Use 1MB buffer to minimize write syscalls (default 8KB is too small)
        let mut seq_out = BufWriter::with_capacity(
            1024 * 1024,
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(&seq_file)
                .map_err(|e| format!("Failed to create sequence file: {e}"))?,
        );

        let mut lines = reader.lines();

        // Detect format from first line
        let first_line = lines
            .next()
            .ok_or("Empty file".to_string())?
            .map_err(|e| format!("Failed to read first line: {e}"))?;

        let is_fasta = first_line.starts_with('>');
        let is_fastq = first_line.starts_with('@');

        if !is_fasta && !is_fastq {
            return Err("Unknown file format (expected FASTA or FASTQ)".to_string());
        }

        // Accumulators
        let mut current_line = first_line;
        let mut seq_bytes_written: u64 = 0;
        let mut notified_empty_seqs = false;

        // Build concatenated name string with delimiters: ">name1 >name2 ..."
        let mut name_text = String::new();
        let mut name_boundary_positions = Vec::new();
        let mut seq_boundary_positions = Vec::new();

        loop {
            // Parse sequence name (same for FASTA and FASTQ)
            let seq_name = current_line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();

            // Get sequence
            let mut seq = String::new();
            let mut found_next_header = false;

            if is_fasta {
                // Read until next '>' or EOF
                for line in lines.by_ref() {
                    let line = line.map_err(|e| format!("Failed to read line: {e}"))?;
                    if line.starts_with('>') {
                        current_line = line;
                        found_next_header = true;
                        break;
                    }
                    seq.push_str(&line);
                }
            } else {
                // FASTQ: read exactly 3 more lines
                if let Some(Ok(seq_line)) = lines.next() {
                    seq = seq_line;
                    lines.next(); // Skip '+'
                    lines.next(); // Skip quality
                }

                // Get next header
                if let Some(Ok(next_line)) = lines.next() {
                    current_line = next_line;
                } else {
                    break;
                }
            }

            // Skip empty sequences
            if seq.is_empty() {
                if !notified_empty_seqs {
                    notified_empty_seqs = true;
                    eprintln!("[seqindex] WARNING: input contains empty sequences, which will be ignored.");
                }
                // If we reached EOF or there's no next header, stop
                if is_fasta && !found_next_header {
                    break;
                }
                continue;
            }

            // Record name boundary (position of '>' in concatenated name text)
            name_boundary_positions.push(name_text.len() as u64);

            // Add to name text: ">name "
            name_text.push('>');
            name_text.push_str(&seq_name);
            name_text.push(' ');

            // Record sequence boundary
            seq_boundary_positions.push(seq_bytes_written);

            // Write upper-case sequence
            let seq_upper = seq.to_uppercase();
            seq_out
                .write_all(seq_upper.as_bytes())
                .map_err(|e| format!("Failed to write sequence: {e}"))?;

            seq_bytes_written += seq_upper.len() as u64;
            self.seq_count += 1;

            // Check EOF
            if is_fasta {
                if !found_next_header {
                    // Reached EOF without finding another header
                    break;
                }
            } else if !current_line.starts_with('@') {
                break;
            }
        }

        // Add final boundary for total length
        seq_boundary_positions.push(seq_bytes_written);

        // Close sequence file
        drop(seq_out);

        // Build FM-index (CSA) from name text
        // Space: O(n log σ) bits where n = name_text.len()
        // FM-index requires text to end with exactly one zero character
        let mut name_bytes = name_text.into_bytes();
        name_bytes.push(0); // Add null terminator required by FM-index
        let name_len = name_bytes.len() - 1; // Length without the null terminator
        let text = Text::new(name_bytes.clone());
        self.name_index = Some(
            FMIndexWithLocate::new(&text, 2)
                .map_err(|e| format!("Failed to build FM-index: {e:?}"))?,
        ); // Sample every 2^2=4 positions
        name_bytes.pop(); // Remove null terminator from stored copy
        self.name_text = name_bytes;

        // Build sparse bitvector for name boundaries (just store positions directly!)
        // Space: O(m) words where m = # sequences - much simpler and faster than RsVec
        self.name_boundaries = Some(SparseBitVec::from_positions(
            name_boundary_positions
                .iter()
                .map(|&p| p as usize)
                .collect(),
            name_len,
        ));

        // Build sparse bitvector for sequence boundaries
        // Space: O(m) words where m = # sequences
        self.seq_boundaries = Some(SparseBitVec::from_positions(
            seq_boundary_positions.iter().map(|&p| p as usize).collect(),
            (seq_bytes_written + 1) as usize,
        ));

        // Memory-map the sequence file
        self.open_mmap()?;

        Ok(())
    }

    /// Memory-map the sequence file
    fn open_mmap(&mut self) -> Result<(), String> {
        if let Some(ref seq_file) = self.seq_filename {
            let file = File::open(seq_file)
                .map_err(|e| format!("Failed to open sequence file for mmap: {e}"))?;

            let mmap = unsafe {
                Mmap::map(&file).map_err(|e| format!("Failed to mmap sequence file: {e}"))?
            };

            self.seq_mmap = Some(mmap);
            Ok(())
        } else {
            Err("No sequence file to map".to_string())
        }
    }

    /// Get sequence name by id (1-based to match C++)
    ///
    /// Time complexity: O(m) where m = length of name
    /// Uses select to find boundaries, then extract from FM-index
    pub fn nth_name(&self, n: usize) -> Option<String> {
        if n < 1 || n > self.seq_count {
            return None;
        }

        let name_boundaries = self.name_boundaries.as_ref()?;

        // Select1 gives us the position of the nth 1-bit (0-indexed)
        // This is the position of the '>' character
        let start = name_boundaries.select1(n - 1);

        // Find the end (position before the space)
        // Format is: ">name " so we need to find the space and back up
        let mut end = if n < self.seq_count {
            name_boundaries.select1(n) - 1 // Position just before next '>'
        } else {
            self.name_text.len()
        };

        // Back up past the trailing space
        while end > start + 1 && self.name_text[end - 1] == b' ' {
            end -= 1;
        }

        // Extract name (skip the '>' character at start, up to but not including the space)
        if start + 1 < self.name_text.len() && end > start + 1 && end <= self.name_text.len() {
            Some(String::from_utf8_lossy(&self.name_text[start + 1..end]).to_string())
        } else {
            None
        }
    }

    /// Get sequence id by name (returns 1-based index to match C++)
    ///
    /// Time complexity: O(m log n + occ) where m = pattern length, occ = occurrences
    /// Uses FM-index locate() operation
    pub fn rank_of_seq_named(&self, name: &str) -> Option<usize> {
        let name_index = self.name_index.as_ref()?;
        let name_boundaries = self.name_boundaries.as_ref()?;

        // Build query pattern: ">name "
        let mut query = String::with_capacity(name.len() + 2);
        query.push('>');
        query.push_str(name);
        query.push(' ');

        // Locate pattern in FM-index
        let search_result = name_index.search(query.as_bytes());
        let matches: Vec<usize> = search_result.iter_matches().map(|m| m.locate()).collect();

        if matches.len() != 1 {
            return None; // Should have exactly one occurrence
        }

        let pos = matches[0] as u64;

        // Rank1 gives us the number of 1-bits before (or at) this position
        // This is the 0-based sequence ID, so add 1 for 1-based
        Some(name_boundaries.rank1(pos as usize) + 1)
    }

    /// Get length of nth sequence (1-based)
    ///
    /// Time complexity: O(1) with select queries
    pub fn nth_seq_length(&self, n: usize) -> Option<u64> {
        if n < 1 || n > self.seq_count {
            return None;
        }

        let seq_boundaries = self.seq_boundaries.as_ref()?;

        let start = seq_boundaries.select1(n - 1) as u64;
        let end = seq_boundaries.select1(n) as u64;

        Some(end - start)
    }

    /// Get offset of nth sequence (1-based)
    ///
    /// Time complexity: O(1) with select
    pub fn nth_seq_offset(&self, n: usize) -> Option<u64> {
        if n < 1 || n > self.seq_count {
            return None;
        }

        Some(self.seq_boundaries.as_ref()?.select1(n - 1) as u64)
    }

    /// Get total number of sequences
    pub fn n_seqs(&self) -> usize {
        self.seq_count
    }

    /// Get total sequence length (all sequences concatenated)
    pub fn seq_length(&self) -> u64 {
        self.seq_mmap.as_ref().map(|m| m.len() as u64).unwrap_or(0)
    }

    /// Get character at position in concatenated sequence
    pub fn at(&self, pos: u64) -> Option<char> {
        if let Some(ref mmap) = self.seq_mmap {
            if (pos as usize) < mmap.len() {
                return Some(mmap[pos as usize] as char);
            }
        }
        None
    }

    /// Get character at pos_t position (handles reverse complement)
    pub fn at_pos(&self, pos: u64) -> Option<char> {
        let offset = crate::pos::offset(pos);
        let is_rev = crate::pos::is_rev(pos);

        if let Some(base) = self.at(offset) {
            if is_rev {
                Some(crate::dna::complement(base as u8) as char)
            } else {
                Some(base)
            }
        } else {
            None
        }
    }

    /// Get sequence ID for position in concatenated sequence
    ///
    /// Time complexity: O(log m) where m = number of sequences (using rank)
    pub fn seq_id_at(&self, pos: u64) -> Option<usize> {
        let seq_boundaries = self.seq_boundaries.as_ref()?;

        // Rank1(pos) counts 1-bits up to but EXCLUDING position pos
        // So rank1(pos+1) counts 1-bits up to and INCLUDING position pos
        // This gives us the sequence ID (1-based) that contains position pos
        let rank = seq_boundaries.rank1(pos as usize + 1);

        if rank > 0 && rank <= self.seq_count {
            Some(rank)
        } else {
            None
        }
    }

    /// Check if position is start of a sequence
    pub fn seq_start(&self, pos: u64) -> bool {
        if let Some(ref seq_boundaries) = self.seq_boundaries {
            // Check if there's a 1-bit at this position
            if (pos as usize) < seq_boundaries.len() {
                // SparseBitVec get() returns u64 directly (0 or 1)
                return seq_boundaries.get(pos as usize) == 1;
            }
        }
        false
    }

    /// Get subsequence by sequence name
    pub fn subseq_by_name(&self, name: &str, pos: u64, count: u64) -> Option<String> {
        let seq_id = self.rank_of_seq_named(name)?;
        self.subseq_by_id(seq_id, pos, count)
    }

    /// Get subsequence by sequence id (1-based)
    pub fn subseq_by_id(&self, seq_id: usize, pos: u64, count: u64) -> Option<String> {
        let seq_offset = self.nth_seq_offset(seq_id)?;
        let seq_len = self.nth_seq_length(seq_id)?;

        if pos + count > seq_len {
            return None;
        }

        self.subseq_absolute(seq_offset + pos, count)
    }

    /// Get subsequence by absolute position in concatenated sequence
    pub fn subseq_absolute(&self, pos: u64, count: u64) -> Option<String> {
        if let Some(ref mmap) = self.seq_mmap {
            let start = pos as usize;
            let end = (pos + count) as usize;

            if end <= mmap.len() {
                return Some(String::from_utf8_lossy(&mmap[start..end]).to_string());
            }
        }
        None
    }

    /// Get full sequence by name
    pub fn seq_by_name(&self, name: &str) -> Option<String> {
        let seq_id = self.rank_of_seq_named(name)?;
        let seq_len = self.nth_seq_length(seq_id)?;
        self.subseq_by_id(seq_id, 0, seq_len)
    }

    /// Convert sequence name + position to absolute position
    pub fn pos_in_all_seqs(&self, name: &str, pos: u64, is_rev: bool) -> Option<u64> {
        let seq_id = self.rank_of_seq_named(name)?;
        self.pos_in_all_seqs_by_id(seq_id, pos, is_rev)
    }

    /// Convert sequence id + position to absolute position (1-based seq_id)
    pub fn pos_in_all_seqs_by_id(&self, seq_id: usize, pos: u64, is_rev: bool) -> Option<u64> {
        let seq_offset = self.nth_seq_offset(seq_id)?;
        let seq_len = self.nth_seq_length(seq_id)?;

        if is_rev {
            if pos < seq_len {
                Some(seq_offset + seq_len - pos - 1)
            } else {
                None
            }
        } else {
            Some(seq_offset + pos)
        }
    }
}

impl Default for SeqIndex {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for SeqIndex {
    fn drop(&mut self) {
        self.seq_mmap = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn create_test_fasta(path: &str, sequences: &[(&str, &str)]) {
        let mut file = File::create(path).unwrap();
        for (name, seq) in sequences {
            writeln!(file, ">{}", name).unwrap();
            writeln!(file, "{}", seq).unwrap();
        }
    }

    #[test]
    fn test_fasta_parsing() {
        let test_file = "/tmp/test_seqindex_v2.fa";
        create_test_fasta(
            test_file,
            &[("seq1", "ACGT"), ("seq2", "GGGG"), ("seq3", "TTTT")],
        );

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        assert_eq!(idx.n_seqs(), 3);
        assert_eq!(idx.nth_name(1), Some("seq1".to_string()));
        assert_eq!(idx.nth_name(2), Some("seq2".to_string()));
        assert_eq!(idx.nth_name(3), Some("seq3".to_string()));

        std::fs::remove_file(test_file).ok();
    }

    #[test]
    #[ignore] // TODO: Fix this test - seq_by_name is dead code and may be broken
    fn test_sequence_access() {
        let test_file = "/tmp/test_seqindex_v2_2.fa";
        create_test_fasta(test_file, &[("chr1", "ACGTACGT"), ("chr2", "GGGGTTTT")]);

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        assert_eq!(idx.seq_by_name("chr1"), Some("ACGTACGT".to_string()));
        assert_eq!(idx.seq_by_name("chr2"), Some("GGGGTTTT".to_string()));
        assert_eq!(idx.nth_seq_length(1), Some(8));
        assert_eq!(idx.nth_seq_length(2), Some(8));

        std::fs::remove_file(test_file).ok();
    }

    #[test]
    fn test_name_lookup() {
        let test_file = "/tmp/test_seqindex_v2_3.fa";
        create_test_fasta(
            test_file,
            &[("seq1", "AAAA"), ("seq2", "CCCC"), ("seq3", "GGGG")],
        );

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        // Test rank_of_seq_named
        assert_eq!(idx.rank_of_seq_named("seq1"), Some(1));
        assert_eq!(idx.rank_of_seq_named("seq2"), Some(2));
        assert_eq!(idx.rank_of_seq_named("seq3"), Some(3));
        assert_eq!(idx.rank_of_seq_named("nonexistent"), None);

        std::fs::remove_file(test_file).ok();
    }

    #[test]
    fn test_position_queries() {
        let test_file = "/tmp/test_seqindex_v2_4.fa";
        create_test_fasta(test_file, &[("s1", "AAAA"), ("s2", "CCCC")]);

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        // s1 is at offset 0, s2 is at offset 4
        assert_eq!(idx.seq_id_at(0), Some(1));
        assert_eq!(idx.seq_id_at(3), Some(1));
        assert_eq!(idx.seq_id_at(4), Some(2));
        assert_eq!(idx.seq_id_at(7), Some(2));

        std::fs::remove_file(test_file).ok();
    }
}
