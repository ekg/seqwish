/// DNA complement lookup table
/// Supports IUPAC nucleotide codes and special characters (GCSA stop/start: $ and #)
const DNA_COMPLEMENT: [u8; 256] = [
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 8
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 16
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 24
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 32
    b'N', b'N', b'N', b'$', b'#', b'N', b'N', b'N', // 40 GCSA stop/start characters
    b'N', b'N', b'N', b'N', b'N', b'-', b'N', b'N', // 48
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 56
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 64
    b'N', b'T', b'V', b'G', b'H', b'N', b'N', b'C', // 72  (A->T, B->V, C->G, D->H)
    b'D', b'N', b'N', b'M', b'N', b'K', b'N', b'N', // 80  (H->D, K->M, M->K)
    b'N', b'Q', b'Y', b'S', b'A', b'A', b'B',
    b'W', // 88  (Q->Q, R->Y, S->W, T->A, U->A, V->B, W->S)
    b'N', b'R', b'N', b'N', b'N', b'N', b'N', b'N', // 96  (Y->R)
    b'N', b't', b'v', b'g', b'h', b'N', b'N', b'c', // 104 (lowercase)
    b'd', b'N', b'N', b'm', b'N', b'k', b'n', b'N', // 112
    b'N', b'q', b'y', b's', b'a', b'a', b'b', b'w', // 120
    b'N', b'r', b'N', b'N', b'N', b'N', b'N', b'N', // 128
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 136
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 144
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 152
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 160
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 168
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 176
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 184
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 192
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 200
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 208
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 216
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 224
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 232
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 240
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 248
    b'N', b'N', b'N', b'N', b'N', b'N', b'N', b'N', // 256
];

/// Get the complement of a single DNA base
#[inline]
pub fn complement(c: u8) -> u8 {
    DNA_COMPLEMENT[c as usize]
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| DNA_COMPLEMENT[c as usize])
        .collect()
}

/// Reverse complement a DNA sequence in place
pub fn reverse_complement_in_place(seq: &mut [u8]) {
    let len = seq.len();
    let swap_size = len / 2;

    // Swap and complement in one pass
    for i in 0..swap_size {
        let j = len - 1 - i;
        let tmp = seq[i];
        seq[i] = DNA_COMPLEMENT[seq[j] as usize];
        seq[j] = DNA_COMPLEMENT[tmp as usize];
    }

    // Handle middle element if odd length
    if len % 2 == 1 {
        seq[swap_size] = DNA_COMPLEMENT[seq[swap_size] as usize];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_basic() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_complement_lowercase() {
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(b't'), b'a');
        assert_eq!(complement(b'c'), b'g');
        assert_eq!(complement(b'g'), b'c');
        assert_eq!(complement(b'n'), b'n');
    }

    #[test]
    fn test_complement_iupac() {
        // IUPAC ambiguity codes
        assert_eq!(complement(b'R'), b'Y'); // R (A/G) -> Y (C/T)
        assert_eq!(complement(b'Y'), b'R'); // Y (C/T) -> R (A/G)
        assert_eq!(complement(b'S'), b'S'); // S (G/C) -> S (G/C)
        assert_eq!(complement(b'W'), b'W'); // W (A/T) -> W (A/T)
        assert_eq!(complement(b'K'), b'M'); // K (G/T) -> M (A/C)
        assert_eq!(complement(b'M'), b'K'); // M (A/C) -> K (G/T)
        assert_eq!(complement(b'B'), b'V'); // B (C/G/T) -> V (A/C/G)
        assert_eq!(complement(b'D'), b'H'); // D (A/G/T) -> H (A/C/T)
        assert_eq!(complement(b'H'), b'D'); // H (A/C/T) -> D (A/G/T)
        assert_eq!(complement(b'V'), b'B'); // V (A/C/G) -> B (C/G/T)
    }

    #[test]
    fn test_complement_special() {
        // GCSA stop/start characters
        assert_eq!(complement(b'$'), b'#');
        assert_eq!(complement(b'#'), b'$');
        // Gap
        assert_eq!(complement(b'-'), b'-');
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"CGCG"), b"CGCG");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement(b"acgt"), b"acgt");
        assert_eq!(reverse_complement(b"aaaa"), b"tttt");
    }

    #[test]
    fn test_reverse_complement_mixed() {
        assert_eq!(reverse_complement(b"AcGt"), b"aCgT");
        assert_eq!(reverse_complement(b"ATCGN"), b"NCGAT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_reverse_complement_single() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"G"), b"C");
    }

    #[test]
    fn test_reverse_complement_in_place_basic() {
        let mut seq = b"ACGT".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"ACGT");

        let mut seq = b"AAAA".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"TTTT");

        let mut seq = b"ATCG".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"CGAT");
    }

    #[test]
    fn test_reverse_complement_in_place_odd_length() {
        let mut seq = b"ATCGN".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"NCGAT");

        let mut seq = b"A".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"T");
    }

    #[test]
    fn test_reverse_complement_in_place_empty() {
        let mut seq = b"".to_vec();
        reverse_complement_in_place(&mut seq);
        assert_eq!(seq, b"");
    }

    #[test]
    fn test_reverse_complement_equivalence() {
        // Test that in-place and allocating versions produce same result
        let test_cases = vec![
            &b"ACGT"[..],
            &b"AAAA"[..],
            &b"ATCGN"[..],
            &b"CGCGCG"[..],
            &b"A"[..],
            &b""[..],
            &b"ATCGATCGATCG"[..],
        ];

        for &test in &test_cases {
            let rc_alloc = reverse_complement(test);
            let mut seq = test.to_vec();
            reverse_complement_in_place(&mut seq);
            assert_eq!(
                rc_alloc,
                seq,
                "Mismatch for input: {:?}",
                std::str::from_utf8(test)
            );
        }
    }

    #[test]
    fn test_reverse_complement_involution() {
        // RC(RC(seq)) = seq
        let test = b"ATCGATCG";
        let rc = reverse_complement(test);
        let rcrc = reverse_complement(&rc);
        assert_eq!(test, &rcrc[..]);
    }
}
