#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * Opaque handle to CIGAR vector
 */
typedef struct CigarHandle CigarHandle;

/**
 * Opaque handle to IITree
 */
typedef struct IITreeHandle IITreeHandle;

/**
 * Opaque handle to a parsed PAF row
 */
typedef struct PafRowHandle PafRowHandle;

/**
 * Opaque handle to SeqIndex
 */
typedef struct SeqIndexHandle SeqIndexHandle;

/**
 * Opaque handle to a parsed SXS alignment
 */
typedef struct SxsHandle SxsHandle;

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

/**
 * Open a file and memory-map it
 * Returns the file size on success, 0 on error
 * The buffer pointer and file descriptor are written to the provided pointers
 */
uintptr_t mmap_open_rust(const char *filename, char **buf_out, int32_t *fd_out);

/**
 * Close a memory-mapped file
 */
void mmap_close_rust(char *buf, int32_t fd, uintptr_t size);

/**
 * Check if a file exists
 */
bool file_exists(const char *filename);

/**
 * Parse a number with optional suffix (k, m, g)
 */
double handy_parameter(const char *value, double default_value);

/**
 * Get milliseconds since Unix epoch
 */
uint64_t time_since_epoch_ms(void);

/**
 * Parse PAF spec string, calling callback for each (filename, weight) pair
 * Callback signature: void callback(void* user_data, const char* filename, uint64_t weight)
 */
void parse_paf_spec(const char *spec,
                    void *user_data,
                    void (*callback)(void*, const char*, uint64_t));

/**
 * Hash function for match parameters
 */
uint64_t match_hash(uint64_t q, uint64_t t, uint64_t l);

/**
 * Determine if a match should be kept based on sparsification factor
 */
bool keep_sparse(uint64_t q, uint64_t t, uint64_t l, float f);

/**
 * Parse a PAF row from a C string line
 * Returns NULL if parsing fails
 */
struct PafRowHandle *paf_row_parse(const char *line);

/**
 * Free a PAF row handle
 */
void paf_row_free(struct PafRowHandle *handle);

char *paf_row_query_sequence_name(const struct PafRowHandle *handle);

char *paf_row_target_sequence_name(const struct PafRowHandle *handle);

uint64_t paf_row_query_sequence_length(const struct PafRowHandle *handle);

uint64_t paf_row_query_start(const struct PafRowHandle *handle);

uint64_t paf_row_query_end(const struct PafRowHandle *handle);

bool paf_row_query_target_same_strand(const struct PafRowHandle *handle);

uint64_t paf_row_target_sequence_length(const struct PafRowHandle *handle);

uint64_t paf_row_target_start(const struct PafRowHandle *handle);

uint64_t paf_row_target_end(const struct PafRowHandle *handle);

uint64_t paf_row_num_matches(const struct PafRowHandle *handle);

uint64_t paf_row_alignment_block_length(const struct PafRowHandle *handle);

uint16_t paf_row_mapping_quality(const struct PafRowHandle *handle);

struct CigarHandle *paf_row_cigar(const struct PafRowHandle *handle);

/**
 * Create a new empty SXS alignment
 */
struct SxsHandle *sxs_new(void);

/**
 * Parse SXS alignment from array of C strings (lines)
 * Returns NULL if parsing fails
 */
struct SxsHandle *sxs_parse_lines(const char *const *lines, uintptr_t num_lines);

/**
 * Free an SXS handle
 */
void sxs_free(struct SxsHandle *handle);

char *sxs_query_sequence_name(const struct SxsHandle *handle);

char *sxs_target_sequence_name(const struct SxsHandle *handle);

uint64_t sxs_query_start(const struct SxsHandle *handle);

uint64_t sxs_query_end(const struct SxsHandle *handle);

uint64_t sxs_target_start(const struct SxsHandle *handle);

uint64_t sxs_target_end(const struct SxsHandle *handle);

uint64_t sxs_num_matches(const struct SxsHandle *handle);

uint16_t sxs_mapping_quality(const struct SxsHandle *handle);

struct CigarHandle *sxs_cigar(const struct SxsHandle *handle);

bool sxs_is_good(const struct SxsHandle *handle);

bool sxs_is_reverse(const struct SxsHandle *handle);

/**
 * Compact nodes by marking boundaries in the graph
 *
 * # Arguments
 * * `seqidx_handle` - Handle to the seqindex
 * * `graph_size` - Size of the graph sequence
 * * `node_iitree_handle` - Handle to the node iitree
 * * `path_iitree_handle` - Handle to the path iitree
 * * `seq_id_bv` - Pointer to bitvector array (will be modified)
 * * `seq_id_bv_size` - Size of the bitvector
 * * `num_threads` - Number of threads to use
 *
 * # Returns
 * 0 on success, 1 on error
 */
int32_t compact_compact_nodes(const struct SeqIndexHandle *seqidx_handle,
                              uintptr_t graph_size,
                              const struct IITreeHandle *node_iitree_handle,
                              const struct IITreeHandle *path_iitree_handle,
                              uint64_t *seq_id_bv,
                              uintptr_t seq_id_bv_size,
                              uintptr_t num_threads);

/**
 * Compute transitive closures for variation graph construction
 *
 * # Arguments
 * * `seqidx_handle` - Handle to the seqindex
 * * `aln_iitree_handle` - Handle to the alignment iitree
 * * `seq_v_file` - Path to output sequence file
 * * `node_iitree_handle` - Handle to the node iitree
 * * `path_iitree_handle` - Handle to the path iitree
 * * `repeat_max` - Maximum repeat count
 * * `min_repeat_dist` - Minimum repeat distance
 * * `transclose_batch_size` - Batch size for transitive closure
 * * `show_progress` - Whether to show progress messages
 * * `num_threads` - Number of threads to use
 *
 * # Returns
 * The length of the graph sequence, or 0 on error
 */
uintptr_t transclosure_compute(const struct SeqIndexHandle *seqidx_handle,
                               const struct IITreeHandle *aln_iitree_handle,
                               const char *seq_v_file,
                               const struct IITreeHandle *node_iitree_handle,
                               const struct IITreeHandle *path_iitree_handle,
                               uint64_t repeat_max,
                               uint64_t min_repeat_dist,
                               uint64_t transclose_batch_size,
                               bool show_progress,
                               uintptr_t num_threads);

char *version_get_version(void);

char *version_get_release(void);

char *version_get_codename(void);

char *version_get_short(void);

void version_free_string(char *s);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
