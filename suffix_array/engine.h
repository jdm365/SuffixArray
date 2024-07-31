#pragma once

#include <stdio.h>
#include <stdint.h>


#define DEFAULT_PAGE_SIZE (uint64_t)4096
#define NUM_BUCKETS 	  (uint32_t)256

typedef struct buffer_c {
	char* 	  buffer;
	uint32_t  buffer_idx;
	uint32_t  buffer_capacity;
} buffer_c;
void init_buffer_c(buffer_c* char_buffer, uint32_t buffer_capacity);
void append_buffer_c(buffer_c* char_buffer, char c);
void free_buffer_c(buffer_c* char_buffer);

typedef struct buffer_u32 {
	uint32_t* __attribute__ ((aligned(64))) buffer;
	uint32_t  buffer_idx;
	uint32_t  buffer_capacity;
} buffer_u32;
void init_buffer_u32(buffer_u32* char_buffer, uint32_t buffer_capacity);
void append_buffer_u32(buffer_u32* char_buffer, uint32_t value);
void free_buffer_u32(buffer_u32* char_buffer);

typedef struct buffer_bit {
	uint8_t*  buffer;
	uint32_t  buffer_bit_idx;
	uint32_t  buffer_capacity;
} buffer_bit;
void init_buffer_bit(buffer_bit* char_buffer, uint32_t buffer_capacity);
void copy_buffer_bit(buffer_bit* src, buffer_bit* dst);
void set_buffer_bit(buffer_bit* char_buffer, uint32_t idx, uint8_t bit);
void append_buffer_bit(buffer_bit* char_buffer, uint8_t bit);
uint8_t get_buffer_bit(buffer_bit* char_buffer, uint32_t idx);
void free_buffer_bit(buffer_bit* char_buffer);


typedef struct SuffixArray {
	uint32_t*   suffix_array;
	buffer_bit* is_quoted_bitflag;
	uint64_t  	global_byte_start_idx;
	uint64_t  	global_byte_end_idx;
	uint32_t  	max_suffix_length;
	uint32_t  	n;
} SuffixArray;

void init_suffix_array(SuffixArray* suffix_array, uint32_t max_suffix_length);
void init_suffix_array_byte_idxs(
		SuffixArray* suffix_array, 
		uint32_t max_suffix_length,
		uint64_t global_byte_start_idx,
		uint64_t global_byte_end_idx,
		uint32_t n
		);
void free_suffix_array(SuffixArray* suffix_array);


void recursive_bucket_sort(
	const char* str,
	uint32_t* suffix_array,
	uint32_t* temp_suffix_array,
	uint32_t overall_size,
	uint32_t current_size,
	int max_depth,
	int current_depth
);

void read_text_into_buffer(
	const char* filename,
	char** buffer,
	uint64_t* buffer_size
);

void construct_truncated_suffix_array_from_csv_partitioned(
	const char* csv_file,
	uint32_t column_idx,
	SuffixArray* suffix_array
);
void construct_truncated_suffix_array(const char* str, SuffixArray* suffix_array);

typedef struct pair_u32 {
	uint32_t first;
	uint32_t second;
} pair_u32;

typedef struct pair_u64 {
	uint64_t first;
	uint64_t second;
} pair_u64;

pair_u32 get_substring_positions(
    const char* str,
	const SuffixArray* suffix_array,
    const char* substring
);

pair_u32 get_substring_positions_file(
    FILE* file,
	const SuffixArray* suffix_array,
    const char* substring
);

uint32_t get_matching_records(
	const char* str,
	const SuffixArray* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records
);

uint32_t get_matching_records_file(
	const char* filename,
	const SuffixArray* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records
);
