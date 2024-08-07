#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


#define DEFAULT_PAGE_SIZE (uint64_t)4096
#define NUM_BUCKETS 	  (uint32_t)256
#define MAX_COUMNS		  4096

inline uint32_t rfc4180_getline(char** lineptr, uint32_t* n, FILE* stream) {
    if (lineptr == NULL || n == NULL || stream == NULL) {
        return UINT32_MAX;
    }

    uint32_t len = 0;
    int c;
    int in_quotes = 0;

    if (*lineptr == NULL) {
        *n = 128;
        *lineptr = (char *)malloc(*n);
        if (*lineptr == NULL) {
            return UINT32_MAX;
        }
    }

    while ((c = fgetc(stream)) != EOF) {
        if (len + 1 >= *n) {
            *n *= 2;
            char *new_lineptr = (char*)realloc(*lineptr, *n);
            if (new_lineptr == NULL) {
                return UINT32_MAX;
            }
            *lineptr = new_lineptr;
        }

        (*lineptr)[len++] = c;

        if (c == '"') {
            in_quotes = !in_quotes;
        } else if (c == '\n' && !in_quotes) {
            break;
        }
    }

    if (ferror(stream)) {
        return UINT32_MAX;
    }

    if (len == 0 && c == EOF) {
        return UINT32_MAX;
    }

    (*lineptr)[len] = '\0';
    return len;
}


inline void parse_csv_header(
		FILE* fh,
		char** columns,
		uint16_t* num_columns
		) {
	char* line;
	uint32_t num_bytes;
	rfc4180_getline(&line, &num_bytes, fh);
	size_t char_idx = 0;
	size_t last_idx = 0;
	size_t col_idx  = 0;

	columns = (char**)malloc(MAX_COUMNS * sizeof(char*));

	while (1) {
		if (line[char_idx] == ',') {
			columns[col_idx] = (char*)malloc((char_idx - last_idx) * sizeof(char));
			strncpy(
					columns[col_idx],
					&line[last_idx],
					char_idx - last_idx
					);

			++col_idx;
			++char_idx;
			last_idx = char_idx;
		} else if (line[char_idx] == '"') {
			++char_idx;
			while (1) {
				if (line[char_idx] == '"') {
					if (line[char_idx + 1] == '"') {
						char_idx += 2;
						break;
					} else {
						++char_idx;
						columns[col_idx] = (char*)malloc((char_idx - last_idx - 2) * sizeof(char));
						strncpy(
								columns[col_idx],
								&line[last_idx],
								char_idx - last_idx
								);

						++col_idx;
						++char_idx;
						last_idx = char_idx;
						break;
					}
				}
				++char_idx;
			}
		} else if ((line[char_idx] == '\n') || (line[char_idx] == '\0')) {
		    columns[col_idx] = (char*)malloc((char_idx - last_idx) * sizeof(char));
			strncpy(
					columns[col_idx],
					&line[last_idx],
					char_idx - last_idx
					);
			// TODO: To lower.

			++col_idx;
			++char_idx;
			last_idx = char_idx;
			break;
		} else {
			++char_idx;
		}
	}

	*num_columns = col_idx;
	columns = (char**)realloc(columns, (size_t)num_columns);
}


typedef struct buffer_c {
	char* 	  buffer;
	uint32_t  buffer_idx;
	uint32_t  buffer_capacity;
} buffer_c;
void init_buffer_c(buffer_c* char_buffer, uint32_t buffer_capacity);
inline void append_buffer_c(buffer_c* buffer, char c) {
	if (buffer->buffer_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (char*)realloc(buffer->buffer, buffer->buffer_capacity);
	}

	buffer->buffer[buffer->buffer_idx++] = c;
}
void free_buffer_c(buffer_c* char_buffer);

typedef struct buffer_u32 {
	uint32_t* __attribute__ ((aligned(64))) buffer;
	uint32_t  buffer_idx;
	uint32_t  buffer_capacity;
} buffer_u32;
void init_buffer_u32(buffer_u32* char_buffer, uint32_t buffer_capacity);
inline void append_buffer_u32(buffer_u32* buffer, uint32_t value) {
	if (buffer->buffer_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (uint32_t*)realloc(
				buffer->buffer, 
				buffer->buffer_capacity * sizeof(uint32_t)
				);
	}

	buffer->buffer[buffer->buffer_idx++] = value;
}
void free_buffer_u32(buffer_u32* char_buffer);

typedef struct buffer_bit {
	uint8_t*  buffer;
	uint32_t  buffer_bit_idx;
	uint32_t  buffer_capacity;
} buffer_bit;
void init_buffer_bit(buffer_bit* char_buffer, uint32_t buffer_capacity);
void copy_buffer_bit(buffer_bit* src, buffer_bit* dst);
void set_buffer_bit(buffer_bit* char_buffer, uint32_t idx, uint8_t bit);
inline void append_buffer_bit(buffer_bit* buffer, uint8_t bit) {
	uint32_t byte_idx = buffer->buffer_bit_idx / 8;
	uint32_t bit_idx  = buffer->buffer_bit_idx % 8;

	if (byte_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (uint8_t*)realloc(buffer->buffer, buffer->buffer_capacity);
	}

	set_buffer_bit(buffer, buffer->buffer_bit_idx, bit);
	++(buffer->buffer_bit_idx);
}
uint8_t get_buffer_bit(buffer_bit* char_buffer, uint32_t idx);
void free_buffer_bit(buffer_bit* char_buffer);
void write_buffer_bit(const buffer_bit* char_buffer, FILE* file);
void read_buffer_bit(buffer_bit* char_buffer, FILE* file);


typedef struct SuffixArray_struct {
	uint32_t*   suffix_array;
	buffer_bit* is_quoted_bitflag;
	uint64_t  	global_byte_start_idx;
	uint64_t  	global_byte_end_idx;
	uint32_t  	max_suffix_length;
	uint32_t  	n;
} SuffixArray_struct;

void init_suffix_array(SuffixArray_struct* suffix_array, uint32_t max_suffix_length);
void init_suffix_array_byte_idxs(
		SuffixArray_struct* suffix_array, 
		uint32_t max_suffix_length,
		uint64_t global_byte_start_idx,
		uint64_t global_byte_end_idx,
		uint32_t n
		);
void free_suffix_array(SuffixArray_struct* suffix_array);
void read_suffix_array(SuffixArray_struct* suffix_array, const char* filename);
void write_suffix_array(
		const SuffixArray_struct* suffix_array, 
		const char* sa_filename,
		const char* is_quoted_filename
		);

typedef struct SuffixArrayFile {
	FILE*		suffix_array_file;
	FILE*		is_quoted_bitflag_file;
	uint64_t  	global_byte_start_idx;
	uint64_t  	global_byte_end_idx;
	uint32_t  	max_suffix_length;
	uint32_t  	n;
} SuffixArrayFile;
void read_into_suffix_array_file(
		SuffixArrayFile* suffix_array_file, 
		const char* sa_filename,
		const char* is_quoted_filename
		);
void free_suffix_array_file(SuffixArrayFile* suffix_array_file);

typedef struct SuffixArrayIndex {
	// TODO: Make union
	SuffixArray_struct* suffix_arrays;
	SuffixArrayFile* suffix_array_files;

	FILE* file_handle;
	buffer_u32 search_col_idxs;
	uint16_t num_partitions;
	uint16_t num_columns;
} SuffixArrayIndex;

void init_suffix_array_index(
		const char* filename,
		SuffixArrayIndex* suffix_array_index,
		uint32_t max_suffix_length,
		char** search_cols,
		uint32_t num_search_cols
		);

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
	SuffixArray_struct* suffix_array
);
void construct_truncated_suffix_array_from_csv_partitioned_mmap(
	const char* csv_file,
	uint32_t column_idx,
	SuffixArray_struct* suffix_array,
	uint16_t num_columns
);
SuffixArrayIndex construct_truncated_suffix_array_from_csv_partitioned_mmap_full(
	const char* csv_file,
	uint32_t max_suffix_length,
	char**   search_cols,
	uint32_t num_search_cols
);

void construct_truncated_suffix_array(const char* str, SuffixArray_struct* suffix_array);

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
	const SuffixArray_struct* suffix_array,
    const char* substring
);

pair_u32 get_substring_positions_file(
    FILE* file,
	const SuffixArray_struct* suffix_array,
    const char* substring
);

pair_u32 get_substring_positions_file_disk(
    FILE* file,
	const SuffixArrayFile* suffix_array,
    const char* substring
);

uint32_t get_matching_records(
	const char* str,
	const SuffixArray_struct* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records
);

void get_matching_records_file(
	const char* filename,
	const SuffixArray_struct* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records,
	uint32_t* num_matches
);

void parse_line(
		const char* line,
		buffer_bit* is_quoted_bitflag,
		buffer_c* text,
		buffer_u32* suffix_array_mapping,
		uint32_t file_pos,
		uint32_t column_idx
		);
