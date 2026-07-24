// AGC input support: reading sequences from an AGC (.agc) archive must yield exactly the same
// sequence index as reading the FASTA the archive was created from. tests/data/agc/small.agc
// was produced with `agc create` from tests/data/agc/small.fa, so ingesting either must give
// identical names, lengths, and bases.
use seqwish::seqindex::SeqIndex;

/// Read a FASTA into (name, sequence) pairs, applying the same "first whitespace token"
/// name rule the index uses.
fn read_fasta(path: &str) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).expect("read fixture FASTA");
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            out.push((name, String::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.push_str(line.trim());
        }
    }
    out
}

#[test]
fn agc_ingest_matches_fasta() {
    let mut fa = SeqIndex::new();
    fa.build_index("tests/data/agc/small.fa")
        .expect("build index from FASTA");

    let mut agc = SeqIndex::new();
    agc.build_index("tests/data/agc/small.agc")
        .expect("build index from AGC");

    assert!(fa.n_seqs() > 0, "fixture should contain sequences");
    assert_eq!(fa.n_seqs(), agc.n_seqs(), "sequence count differs");

    // Same names and lengths, in the same order.
    for i in 1..=fa.n_seqs() {
        assert_eq!(
            fa.nth_name(i),
            agc.nth_name(i),
            "sequence name differs at index {i}"
        );
        assert_eq!(
            fa.nth_seq_length(i),
            agc.nth_seq_length(i),
            "sequence length differs at index {i}"
        );
    }

    // Same concatenated bases, position by position.
    assert_eq!(
        fa.seq_length(),
        agc.seq_length(),
        "total concatenated length differs"
    );
    for p in 0..fa.seq_length() {
        assert_eq!(
            fa.at(p),
            agc.at(p),
            "base differs at concatenated position {p}"
        );
    }
}

/// AGC identifies a sequence by (sample, contig) while seqwish identifies it by name alone.
/// tests/data/agc/dup.agc holds two samples (dupA, dupB) that BOTH contain chr1 and chr2, the
/// canonical AGC layout (`agc create sampleA.fa sampleB.fa`). Such records must stay distinct:
/// each gets a "contig@sample" name and must carry its own sample's bases, not the bases of
/// whichever sample happens to be listed first.
#[test]
fn agc_duplicate_contig_names_stay_distinct() {
    let mut idx = SeqIndex::new();
    idx.build_index("tests/data/agc/dup.agc")
        .expect("build index from AGC with duplicate contig names");

    // Expected content, taken from the FASTAs the archive was built from.
    let mut expected: Vec<(String, String)> = Vec::new();
    for (sample, path) in [
        ("dupA", "tests/data/agc/dupA.fa"),
        ("dupB", "tests/data/agc/dupB.fa"),
    ] {
        for (name, seq) in read_fasta(path) {
            expected.push((format!("{name}@{sample}"), seq));
        }
    }

    assert_eq!(idx.n_seqs(), expected.len(), "sequence count differs");

    for (name, seq) in &expected {
        let id = idx
            .rank_of_seq_named(name)
            .unwrap_or_else(|| panic!("sequence '{name}' not found (names must be unique)"));
        let len = idx.nth_seq_length(id).expect("length of sequence");
        assert_eq!(len as usize, seq.len(), "length differs for '{name}'");
        let got = idx.subseq_by_id(id, 0, len).expect("bases of sequence");
        assert_eq!(&got, seq, "bases differ for '{name}' (wrong sample?)");
    }
}
