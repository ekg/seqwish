#!/usr/bin/env bash
# Test the Rust seqwish implementation against all HLA datasets

set -e

SEQWISH_RUST="./seqwish-rs/target/release/seqwish"
TEMP_DIR="/tmp/seqwish-test-$$"
RESULTS_DIR="/tmp/seqwish-results-$$"

mkdir -p "$TEMP_DIR" "$RESULTS_DIR"

# Test datasets
DATASETS=(
    "V-352962"
    "E-3133"
    "J-3137"
    "K-3138"
    "H-3136"
    "L-3139"
    "F-3134"
    "G-3135"
    "A-3105"
    "B-3106"
    "C-3107"
)

echo "Testing Rust seqwish implementation"
echo "===================================="
echo ""

PASSED=0
FAILED=0

for dataset in "${DATASETS[@]}"; do
    echo -n "Testing $dataset... "

    SEQ="test/HLA/${dataset}.fa.gz"
    PAF="test/HLA/${dataset}.paf.gz"
    OUT="$RESULTS_DIR/${dataset}.gfa"
    EXPECTED_MD5=$(cat "test/HLA/${dataset}.fa.gz.gfa.md5")

    # Run seqwish
    if $SEQWISH_RUST -s "$SEQ" -p "$PAF" -g "$OUT" -b "$TEMP_DIR" > /dev/null 2>&1; then
        # Check MD5
        ACTUAL_MD5=$(md5sum "$OUT" | cut -f1 -d' ')

        if [ "$ACTUAL_MD5" = "$EXPECTED_MD5" ]; then
            echo "✅ PASS"
            ((PASSED++))
        else
            echo "❌ FAIL (checksum mismatch)"
            echo "  Expected: $EXPECTED_MD5"
            echo "  Got:      $ACTUAL_MD5"
            ((FAILED++))
        fi
    else
        echo "❌ FAIL (execution error)"
        ((FAILED++))
    fi
done

echo ""
echo "===================================="
echo "Results: $PASSED passed, $FAILED failed"
echo "===================================="

# Cleanup
rm -rf "$TEMP_DIR" "$RESULTS_DIR"

if [ $FAILED -eq 0 ]; then
    exit 0
else
    exit 1
fi
