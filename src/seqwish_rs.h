#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * Opaque handle to CIGAR vector
 */
typedef struct CigarHandle CigarHandle;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Returns the version string of the Rust component
 */
const char *seqwish_rust_version(void);

/**
 * Simple test function to verify FFI is working
 */
int32_t seqwish_rust_add(int32_t a, int32_t b);

/**
 * Create a temporary file. Returns a C string that must be freed with temp_file_free_string.
 * Returns NULL on error.
 */
char *temp_file_create(const char *base, const char *suffix);

/**
 * Remove a temporary file
 */
void temp_file_remove(const char *filename);

/**
 * Set temp directory
 */
void temp_file_set_dir(const char *dir);

/**
 * Get temp directory. Returns a C string that must be freed with temp_file_free_string.
 */
char *temp_file_get_dir(void);

/**
 * Set whether to keep temp files
 */
void temp_file_set_keep_temp(bool setting);

/**
 * Free a string returned by temp_file functions
 */
void temp_file_free_string(char *s);

/**
 * Create a position from offset and orientation
 */
uint64_t pos_make_pos_t(uint64_t offset, bool is_rev);

/**
 * Extract offset from position
 */
uint64_t pos_offset(uint64_t pos);

/**
 * Check if position is reverse
 */
bool pos_is_rev(uint64_t pos);

/**
 * Increment position
 */
void pos_incr_pos(uint64_t *pos);

/**
 * Increment position by N
 */
void pos_incr_pos_by(uint64_t *pos, uintptr_t by);

/**
 * Decrement position
 */
void pos_decr_pos(uint64_t *pos);

/**
 * Decrement position by N
 */
void pos_decr_pos_by(uint64_t *pos, uintptr_t by);

/**
 * Reverse position orientation
 */
uint64_t pos_rev_pos_t(uint64_t pos);

/**
 * Convert position to string (returns C string that must be freed)
 */
char *pos_to_string_c(uint64_t pos);

/**
 * Get complement of a single DNA base
 */
uint8_t dna_complement(uint8_t c);

/**
 * Reverse complement a DNA sequence (allocates new string that must be freed)
 */
void dna_reverse_complement(const char *seq, uintptr_t len, char *out);

/**
 * Reverse complement a DNA sequence in place
 */
void dna_reverse_complement_in_place(char *seq, uintptr_t len);

/**
 * Parse CIGAR string and return handle to CIGAR vector
 * Returns NULL on error. Must be freed with cigar_free.
 */
struct CigarHandle *cigar_from_string(const char *s);

/**
 * Convert CIGAR vector to string
 * Returns C string that must be freed with temp_file_free_string
 */
char *cigar_to_string(const struct CigarHandle *handle);

/**
 * Get number of operations in CIGAR
 */
uintptr_t cigar_length(const struct CigarHandle *handle);

/**
 * Get operation at index
 * Returns false if index out of bounds
 */
bool cigar_get_op(const struct CigarHandle *handle,
                  uintptr_t index,
                  uint64_t *len_out,
                  uint8_t *op_out);

/**
 * Free CIGAR handle
 */
void cigar_free(struct CigarHandle *handle);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
