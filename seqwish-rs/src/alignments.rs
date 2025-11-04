/// Hash function for match parameters
/// Uses a Wang hash-like mixing algorithm
pub fn match_hash(q: u64, t: u64, l: u64) -> u64 {
    let mut seed = q | t | l;
    seed ^= q.wrapping_add(0x9e3779b97f4a7c15)
        .wrapping_add(seed << 17)
        .wrapping_add(seed >> 9);
    seed ^= t.wrapping_add(0x9e3779b97f4a7c15)
        .wrapping_add(seed << 7)
        .wrapping_add(seed >> 23);
    seed ^= l.wrapping_add(0x9e3779b97f4a7c15)
        .wrapping_add(seed << 9)
        .wrapping_add(seed >> 2);
    seed
}

/// Determine if a match should be kept based on sparsification factor
/// Returns true if the match hash is below the threshold determined by f
pub fn keep_sparse(q: u64, t: u64, l: u64, f: f32) -> bool {
    // Match C++ behavior: compare hash as f64 against max * f
    (match_hash(q, t, l) as f64) < (u64::MAX as f64 * f as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_match_hash_deterministic() {
        // Same inputs should give same hash
        let hash1 = match_hash(100, 200, 50);
        let hash2 = match_hash(100, 200, 50);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_match_hash_different_inputs() {
        // Different inputs should give different hashes (usually)
        let hash1 = match_hash(100, 200, 50);
        let hash2 = match_hash(100, 200, 51);
        assert_ne!(hash1, hash2);

        let hash3 = match_hash(100, 201, 50);
        assert_ne!(hash1, hash3);

        let hash4 = match_hash(101, 200, 50);
        assert_ne!(hash1, hash4);
    }

    #[test]
    fn test_match_hash_zero_inputs() {
        let hash = match_hash(0, 0, 0);
        // Verify it matches the C++ implementation
        assert_eq!(hash, 668627425924508912);
    }

    #[test]
    fn test_keep_sparse_always_keep() {
        // f = 1.0 should keep everything
        assert!(keep_sparse(100, 200, 50, 1.0));
        assert!(keep_sparse(0, 0, 0, 1.0));
        assert!(keep_sparse(u64::MAX, u64::MAX, u64::MAX, 1.0));
    }

    #[test]
    fn test_keep_sparse_never_keep() {
        // f = 0.0 should keep nothing
        assert!(!keep_sparse(100, 200, 50, 0.0));
        assert!(!keep_sparse(0, 0, 0, 0.0));
        assert!(!keep_sparse(1, 1, 1, 0.0));
    }

    #[test]
    fn test_keep_sparse_threshold_behavior() {
        // The hash function produces values in a clustered range
        // Test that the threshold comparison works correctly
        let hash_val = match_hash(100, 200, 50);

        // With f = 1.0, everything should pass
        assert!((hash_val as f64) < (u64::MAX as f64 * 1.0));

        // With f = 0.0, nothing should pass
        assert!(!((hash_val as f64) < (u64::MAX as f64 * 0.0)));

        // With intermediate values, it depends on the hash
        // Just verify the comparison logic works
        let small_f = 0.0001;
        let large_f = 0.9999;
        let _small_result = keep_sparse(100, 200, 50, small_f);
        let large_result = keep_sparse(100, 200, 50, large_f);
        assert!(large_result); // Large f should keep most things
    }

    #[test]
    fn test_keep_sparse_consistency() {
        // Same inputs should give same result
        let result1 = keep_sparse(12345, 67890, 111, 0.3);
        let result2 = keep_sparse(12345, 67890, 111, 0.3);
        assert_eq!(result1, result2);
    }
}
