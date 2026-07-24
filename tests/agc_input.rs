// AGC input support: reading sequences from an AGC (.agc) archive must yield exactly the same
// sequence index as reading the FASTA the archive was created from. tests/data/agc/small.agc
// was produced with `agc create` from tests/data/agc/small.fa, so ingesting either must give
// identical names, lengths, and bases.
use seqwish::seqindex::SeqIndex;

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
