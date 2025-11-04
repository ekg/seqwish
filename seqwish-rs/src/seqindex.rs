use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use memmap2::Mmap;
use flate2::read::GzDecoder;

/// Sequence index for FASTA/FASTQ files
/// Provides random access to sequences by name or position
pub struct SeqIndex {
    /// Base filename for temp files
    base_filename: Option<PathBuf>,

    /// Path to concatenated sequence file
    seq_filename: Option<PathBuf>,

    /// Memory-mapped sequence data
    seq_mmap: Option<Mmap>,

    /// Sequence name -> id mapping
    name_to_id: HashMap<String, usize>,

    /// Sequence id -> name mapping
    id_to_name: Vec<String>,

    /// Sequence start offsets in concatenated file
    seq_offsets: Vec<u64>,

    /// Total number of sequences
    seq_count: usize,
}

impl SeqIndex {
    /// Create a new empty sequence index
    pub fn new() -> Self {
        SeqIndex {
            base_filename: None,
            seq_filename: None,
            seq_mmap: None,
            name_to_id: HashMap::new(),
            id_to_name: Vec::new(),
            seq_offsets: Vec::new(),
            seq_count: 0,
        }
    }

    /// Build index from FASTA or FASTQ file
    pub fn build_index(&mut self, filename: &str) -> Result<(), String> {
        // Create temp files
        let seq_file = std::env::temp_dir().join(format!("seqwish-{}.sqq", std::process::id()));
        self.seq_filename = Some(seq_file.clone());

        // Open input file (with optional gzip support)
        let file = File::open(filename)
            .map_err(|e| format!("Failed to open {}: {}", filename, e))?;

        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        // Open output file for sequences
        let mut seq_out = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(&seq_file)
            .map_err(|e| format!("Failed to create sequence file: {}", e))?;

        let mut lines = reader.lines();

        // Detect format from first line
        let first_line = lines.next()
            .ok_or("Empty file".to_string())?
            .map_err(|e| format!("Failed to read first line: {}", e))?;

        let is_fasta = first_line.starts_with('>');
        let is_fastq = first_line.starts_with('@');

        if !is_fasta && !is_fastq {
            return Err(format!("Unknown file format (expected FASTA or FASTQ)"));
        }

        let mut current_line = first_line;
        let mut seq_bytes_written: u64 = 0;
        let mut notified_empty_seqs = false;

        loop {
            // Parse sequence name (first_line is the header)
            let seq_name = if is_fasta {
                current_line[1..].split_whitespace().next()
                    .unwrap_or("").to_string()
            } else {
                current_line[1..].split_whitespace().next()
                    .unwrap_or("").to_string()
            };

            // Get sequence
            let mut seq = String::new();

            if is_fasta {
                // Read until next '>' or EOF
                for line in lines.by_ref() {
                    let line = line.map_err(|e| format!("Failed to read line: {}", e))?;
                    if line.starts_with('>') {
                        current_line = line;
                        break;
                    }
                    seq.push_str(&line);
                }
            } else {
                // FASTQ: read exactly 3 more lines (seq, +, qual)
                if let Some(Ok(seq_line)) = lines.next() {
                    seq = seq_line;
                    lines.next(); // Skip '+' line
                    lines.next(); // Skip quality line
                }

                // Get next header
                if let Some(Ok(next_line)) = lines.next() {
                    current_line = next_line;
                } else {
                    // EOF
                    break;
                }
            }

            // Skip empty sequences
            if seq.is_empty() {
                if !notified_empty_seqs {
                    notified_empty_seqs = true;
                    eprintln!("[seqindex] WARNING: input contains empty sequences, which will be ignored.");
                }

                if is_fasta && lines.by_ref().next().is_none() {
                    break;
                }
                continue;
            }

            // Store sequence name and offset
            let seq_id = self.seq_count;
            self.name_to_id.insert(seq_name.clone(), seq_id);
            self.id_to_name.push(seq_name);
            self.seq_offsets.push(seq_bytes_written);

            // Write upper-case sequence
            let seq_upper = seq.to_uppercase();
            seq_out.write_all(seq_upper.as_bytes())
                .map_err(|e| format!("Failed to write sequence: {}", e))?;

            seq_bytes_written += seq_upper.len() as u64;
            self.seq_count += 1;

            // Check if we've reached EOF
            if is_fasta {
                if !current_line.starts_with('>') {
                    break;
                }
            } else {
                if !current_line.starts_with('@') {
                    break;
                }
            }
        }

        // Add final offset for total length
        self.seq_offsets.push(seq_bytes_written);

        // Close and memory-map the sequence file
        drop(seq_out);
        self.open_mmap()?;

        Ok(())
    }

    /// Memory-map the sequence file
    fn open_mmap(&mut self) -> Result<(), String> {
        if let Some(ref seq_file) = self.seq_filename {
            let file = File::open(seq_file)
                .map_err(|e| format!("Failed to open sequence file for mmap: {}", e))?;

            let mmap = unsafe {
                Mmap::map(&file)
                    .map_err(|e| format!("Failed to mmap sequence file: {}", e))?
            };

            self.seq_mmap = Some(mmap);
            Ok(())
        } else {
            Err("No sequence file to map".to_string())
        }
    }

    /// Get sequence name by id
    pub fn nth_name(&self, n: usize) -> Option<&str> {
        self.id_to_name.get(n).map(|s| s.as_str())
    }

    /// Get sequence id by name
    pub fn rank_of_seq_named(&self, name: &str) -> Option<usize> {
        self.name_to_id.get(name).copied()
    }

    /// Get length of nth sequence
    pub fn nth_seq_length(&self, n: usize) -> Option<u64> {
        if n < self.seq_count {
            Some(self.seq_offsets[n + 1] - self.seq_offsets[n])
        } else {
            None
        }
    }

    /// Get offset of nth sequence
    pub fn nth_seq_offset(&self, n: usize) -> Option<u64> {
        self.seq_offsets.get(n).copied()
    }

    /// Get total number of sequences
    pub fn n_seqs(&self) -> usize {
        self.seq_count
    }

    /// Get total sequence length (all sequences concatenated)
    pub fn seq_length(&self) -> u64 {
        self.seq_offsets.last().copied().unwrap_or(0)
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
        // Extract offset and reverse flag from pos_t
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
    pub fn seq_id_at(&self, pos: u64) -> Option<usize> {
        // Binary search for the sequence containing this position
        match self.seq_offsets.binary_search(&pos) {
            Ok(idx) => Some(idx),
            Err(idx) => {
                if idx > 0 && idx <= self.seq_count {
                    Some(idx - 1)
                } else {
                    None
                }
            }
        }
    }

    /// Check if position is start of a sequence
    pub fn seq_start(&self, pos: u64) -> bool {
        self.seq_offsets.iter().any(|&offset| offset == pos)
    }

    /// Get subsequence by sequence name
    pub fn subseq_by_name(&self, name: &str, pos: u64, count: u64) -> Option<String> {
        let seq_id = self.rank_of_seq_named(name)?;
        self.subseq_by_id(seq_id, pos, count)
    }

    /// Get subsequence by sequence id
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
                return Some(
                    String::from_utf8_lossy(&mmap[start..end]).to_string()
                );
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

    /// Convert sequence id + position to absolute position
    pub fn pos_in_all_seqs_by_id(&self, seq_id: usize, pos: u64, is_rev: bool) -> Option<u64> {
        let seq_offset = self.nth_seq_offset(seq_id)?;
        let seq_len = self.nth_seq_length(seq_id)?;

        if is_rev {
            // Reverse position
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
        // Clean up memory map
        self.seq_mmap = None;

        // Optionally remove temp files
        // (C++ version has temp_file::keep_temp flag we could respect)
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
        let test_file = "/tmp/test_seqindex.fa";
        create_test_fasta(test_file, &[
            ("seq1", "ACGT"),
            ("seq2", "GGGG"),
            ("seq3", "TTTT"),
        ]);

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        assert_eq!(idx.n_seqs(), 3);
        assert_eq!(idx.nth_name(0), Some("seq1"));
        assert_eq!(idx.nth_name(1), Some("seq2"));
        assert_eq!(idx.nth_name(2), Some("seq3"));

        std::fs::remove_file(test_file).ok();
    }

    #[test]
    fn test_sequence_access() {
        let test_file = "/tmp/test_seqindex2.fa";
        create_test_fasta(test_file, &[
            ("chr1", "ACGTACGT"),
            ("chr2", "GGGGTTTT"),
        ]);

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        assert_eq!(idx.seq_by_name("chr1"), Some("ACGTACGT".to_string()));
        assert_eq!(idx.seq_by_name("chr2"), Some("GGGGTTTT".to_string()));
        assert_eq!(idx.nth_seq_length(0), Some(8));
        assert_eq!(idx.nth_seq_length(1), Some(8));

        std::fs::remove_file(test_file).ok();
    }

    #[test]
    fn test_position_queries() {
        let test_file = "/tmp/test_seqindex3.fa";
        create_test_fasta(test_file, &[
            ("s1", "AAAA"),
            ("s2", "CCCC"),
        ]);

        let mut idx = SeqIndex::new();
        idx.build_index(test_file).unwrap();

        // s1 is at offset 0, s2 is at offset 4
        assert_eq!(idx.seq_id_at(0), Some(0));
        assert_eq!(idx.seq_id_at(3), Some(0));
        assert_eq!(idx.seq_id_at(4), Some(1));
        assert_eq!(idx.seq_id_at(7), Some(1));

        std::fs::remove_file(test_file).ok();
    }
}
