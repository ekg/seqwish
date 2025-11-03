#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

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

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
