#pragma once

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
	std::vector<uint32_t>& suffix_array,
	uint32_t n,
	uint32_t max_suffix_length
);

void construct_truncated_suffix_array_from_csv_partitioned(
	const char* csv_file,
	uint32_t column_idx,
	std::vector<uint32_t>& suffix_array,
	uint32_t* suffix_array_size,
	uint32_t max_suffix_length,
	uint64_t start_idx,
	uint64_t& end_idx
);

std::pair<uint32_t, uint32_t> get_substring_positions(
    const char* str,
    uint32_t* suffix_array,
    uint64_t n,
    const char* substring
);

std::pair<uint64_t, uint64_t> get_substring_positions_file(
    FILE* file,
	uint64_t byte_offset,
    uint32_t* suffix_array,
    uint32_t n,
    const char* substring
);

std::vector<std::string> get_matching_records(
	const char* str,
	uint32_t* suffix_array,
	uint32_t n,
	const char* substring,
	int k 
);

std::vector<std::string> get_matching_records_file(
	const char* filename,
	uint32_t* suffix_array,
	uint32_t n,
	uint64_t byte_offset,
	const char* substring,
	int k 
);
