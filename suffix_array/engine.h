#include <stdint.h>
#include <vector>

void recursive_bucket_sort(
	const char* str,
	uint32_t* suffix_array,
	uint32_t* temp_suffix_array,
	int string_length,
	uint64_t n,
	int max_depth,
	int current_depth
);

void read_text_into_buffer(
	const char* filename,
	char** buffer,
	uint64_t* buffer_size
);

void construct_truncated_suffix_array(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	uint32_t max_suffix_length
);

std::vector<uint32_t> get_substring_positions(
    const char* str,
    uint32_t* suffix_array,
    uint64_t n,
    const char* substring
);

std::vector<uint32_t> get_matching_indices(
	const char* str,
	uint32_t* suffix_array,
	uint32_t* suffix_array_idxs,
	uint64_t n,
	const char* substring,
	int k 
);
