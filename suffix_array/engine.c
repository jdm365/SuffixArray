#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <immintrin.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>

#include "engine.h"


#ifdef NDEBUG
# define NDEBUG_DISABLED
# undef  NDEBUG
#endif
#include <assert.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

inline void swap(uint32_t a, uint32_t b) {
	a ^= b;
	b ^= a;
	a ^= b;
}

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


inline void parse_line(
		const char* line,
		buffer_bit* is_quoted_bitflag,
		buffer_c* text,
		buffer_u32* suffix_array_mapping,
		uint32_t file_pos,
		uint32_t column_idx
		) {

	uint32_t char_idx = 0;
	uint32_t col_idx  = 0;

	uint32_t line_length = strlen(line);

	uint16_t _search_col_idx = 0;

	// Iterate of line chars until we get to relevant column.
	while (col_idx != column_idx) {
		if (line[char_idx] == '\n') {
			printf("Newline found before end.\n");
			printf("Search col idx: %d\n", column_idx);
			printf("Col idx: %d\n", col_idx);
			printf("Char idx: %d\n", char_idx);
			printf("Line: %s", line);
			exit(1);
		}

		if (line[char_idx] == '"') {
			// Skip to next unescaped quote
			++char_idx;

			while (1) {
				if (line[char_idx] == '"') {
					if (line[char_idx + 1] == '"') {
						char_idx += 2;
						continue;
					} 
					else {
						++char_idx;
						break;
					}
				}
				++char_idx;
			}
		}

		if (line[char_idx] == ',') ++col_idx;
		++char_idx;
	}

	uint8_t quoted_field = (line[char_idx] == '"');

	if (quoted_field) {
		++char_idx;

		while (1) {
			if (char_idx > 1048576) {
				printf("Error: Column not found.\n");
				exit(1);
			}

			if (line[char_idx] == '"' && line[char_idx + 1] == '"') {
				append_buffer_c(text, line[char_idx]);
				append_buffer_u32(suffix_array_mapping, file_pos + char_idx);
				append_buffer_bit(is_quoted_bitflag, 1);
				char_idx += 2;
				continue;
			}

			if (
					line[char_idx] == '"' 
						&& 
					(
						line[char_idx + 1] == '\n' 
							|| 
						line[char_idx + 1] == ',' 
							|| 
						line[char_idx + 1] == '\0' 
					)
				) {
				++char_idx;
				break;
			}

			append_buffer_c(text, tolower(line[char_idx]));
			append_buffer_u32(suffix_array_mapping, file_pos + char_idx);
			append_buffer_bit(is_quoted_bitflag, 1);
			++char_idx;
		}
	} else {

		while (1) {
			if (char_idx > 1048576) {
				printf("Error: Column not found.\n");
				exit(1);
			}

			if (line[char_idx] == '\n' || line[char_idx] == '\0' || line[char_idx] == ',') {
				break;
			}

			append_buffer_c(text, tolower(line[char_idx]));
			append_buffer_u32(suffix_array_mapping, file_pos + char_idx);
			append_buffer_bit(is_quoted_bitflag, 0);
			++char_idx;
		}
	}

	append_buffer_c(text, '\n');
	append_buffer_u32(suffix_array_mapping, file_pos + char_idx);
	append_buffer_bit(is_quoted_bitflag, quoted_field);
}

void init_buffer_c(buffer_c* buffer, uint32_t buffer_capacity) {
	if (buffer == NULL) {
		printf("Error: Buffer is NULL.\n");
		exit(1);
	}
	buffer->buffer_capacity = buffer_capacity;
	buffer->buffer_idx = 0;
	buffer->buffer = (char*)malloc(buffer_capacity);
}

void append_buffer_c(buffer_c* buffer, char c) {
	if (buffer->buffer_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (char*)realloc(buffer->buffer, buffer->buffer_capacity);
	}

	buffer->buffer[buffer->buffer_idx++] = c;
}

void free_buffer_c(buffer_c* buffer) {
	free(buffer->buffer);
}


void init_buffer_u32(buffer_u32* buffer, uint32_t buffer_capacity) {
	if (buffer == NULL) {
		printf("Error: Buffer is NULL.\n");
		exit(1);
	}
	buffer->buffer_capacity = buffer_capacity;
	buffer->buffer_idx = 0;
	buffer->buffer = (uint32_t*)malloc(buffer_capacity * sizeof(uint32_t));
}

void append_buffer_u32(buffer_u32* buffer, uint32_t value) {
	if (buffer->buffer_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (uint32_t*)realloc(
				buffer->buffer, 
				buffer->buffer_capacity * sizeof(uint32_t)
				);
	}

	buffer->buffer[buffer->buffer_idx++] = value;
}

void free_buffer_u32(buffer_u32* buffer) {
	free(buffer->buffer);
}


void init_buffer_bit(buffer_bit* buffer, uint32_t buffer_capacity) {
	if (buffer == NULL) {
		printf("Error: Buffer is NULL.\n");
		exit(1);
	}
	buffer->buffer_capacity = buffer_capacity;
	buffer->buffer_bit_idx = 0;
	buffer->buffer = (uint8_t*)malloc(buffer_capacity);
	memset(buffer->buffer, 0, buffer_capacity);
}

void copy_buffer_bit(buffer_bit* src, buffer_bit* dst) {
	memcpy(dst->buffer, src->buffer, src->buffer_capacity);
	dst->buffer_capacity = src->buffer_capacity;
	dst->buffer_bit_idx  = src->buffer_bit_idx;
}

uint8_t get_buffer_bit(buffer_bit* buffer, uint32_t idx) {
	uint32_t byte_idx = idx / 8;
	uint32_t bit_idx  = idx % 8;

	return (buffer->buffer[byte_idx] >> bit_idx) & 1;
}

void set_buffer_bit(buffer_bit* buffer, uint32_t idx, uint8_t bit) {
	uint32_t byte_idx = idx / 8;
	uint32_t bit_idx  = idx % 8;

	if (byte_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (uint8_t*)realloc(buffer->buffer, buffer->buffer_capacity);
	}

	// buffer->buffer[byte_idx] |= (bit << bit_idx);
	if (bit) {
		buffer->buffer[byte_idx] |= (1 << bit_idx);
	} else {
		buffer->buffer[byte_idx] &= ~(1 << bit_idx);
	}
}

static inline void print_bits(uint8_t byte) {
	for (uint32_t i = 0; i < 8; ++i) {
		printf("%d ", (byte >> i) & 1);
	}
}

static inline void print_text_near_idx(
		FILE* file,
		uint64_t idx
		) {
	char line[64];
	fseek(file, idx - 32, SEEK_SET);
	fread(line, 1, 64, file);
	printf(" Line: %.*s", 64, line);
}

void append_buffer_bit(buffer_bit* buffer, uint8_t bit) {
	uint32_t byte_idx = buffer->buffer_bit_idx / 8;
	uint32_t bit_idx  = buffer->buffer_bit_idx % 8;

	if (byte_idx == buffer->buffer_capacity) {
		buffer->buffer_capacity *= 2;
		buffer->buffer = (uint8_t*)realloc(buffer->buffer, buffer->buffer_capacity);
	}

	set_buffer_bit(buffer, buffer->buffer_bit_idx, bit);
	++(buffer->buffer_bit_idx);
}



void free_buffer_bit(buffer_bit* buffer) {
	free(buffer->buffer);
}

void init_suffix_array(
		SuffixArray* suffix_array, 
		uint32_t max_suffix_length
		) {
	if (suffix_array == NULL) {
		printf("Error: Suffix array is NULL.\n");
		exit(1);
	}
	suffix_array->max_suffix_length = max_suffix_length;
	suffix_array->n = 0;
	suffix_array->suffix_array = NULL;
	suffix_array->global_byte_start_idx = 0;
	suffix_array->global_byte_end_idx   = 0;
	suffix_array->is_quoted_bitflag = (buffer_bit*)malloc(sizeof(buffer_bit));
	init_buffer_bit(suffix_array->is_quoted_bitflag, 1024);
}

void init_suffix_array_byte_idxs(
		SuffixArray* suffix_array, 
		uint32_t max_suffix_length,
		uint64_t global_byte_start_idx,
		uint64_t global_byte_end_idx,
		uint32_t n
		) {
	if (suffix_array == NULL) {
		printf("Error: Suffix array is NULL.\n");
		exit(1);
	}
	suffix_array->max_suffix_length = max_suffix_length;
	suffix_array->n = n;
	suffix_array->suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));
	suffix_array->is_quoted_bitflag = (buffer_bit*)malloc(sizeof(buffer_bit));
	init_buffer_bit(suffix_array->is_quoted_bitflag, n);
	suffix_array->global_byte_start_idx = global_byte_start_idx;
	suffix_array->global_byte_end_idx   = global_byte_end_idx;
}

void free_suffix_array(SuffixArray* suffix_array) {
	free(suffix_array->suffix_array);
	free_buffer_bit(suffix_array->is_quoted_bitflag);
}

void construct_truncated_suffix_array_from_csv_partitioned(
	const char* csv_file,
	uint32_t column_idx,
	SuffixArray* suffix_array
) {

	const uint32_t TWO_GB = (uint32_t)2 * (uint32_t)1024 * (uint32_t)1024 * (uint32_t)1024;

	// TODO: Adjust for multithreading correctness.
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	// Read and parse the CSV file.
	FILE* file = fopen(csv_file, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	char*    line 		= NULL;
	uint32_t len  		= 0;
	uint32_t file_pos 	= suffix_array->global_byte_start_idx;
	uint32_t bytes_read = 0;
	uint32_t read;

	uint32_t num_lines_guess = TWO_GB / (suffix_array->max_suffix_length * 8);
	uint32_t num_bytes_guess = num_lines_guess * suffix_array->max_suffix_length;
	fseek(file, suffix_array->global_byte_start_idx, SEEK_SET);

	buffer_c text;
	init_buffer_c(&text, num_bytes_guess);

	buffer_u32 suffix_array_mapping;
	init_buffer_u32(&suffix_array_mapping, num_bytes_guess);

	init_buffer_bit(suffix_array->is_quoted_bitflag, num_lines_guess);

	while (bytes_read < TWO_GB) {
		read = rfc4180_getline(&line, &len, file);
		if (read == UINT32_MAX) break;

		parse_line(
			line,
			suffix_array->is_quoted_bitflag,
			&text,
			&suffix_array_mapping,
			bytes_read,
			column_idx
		);

		bytes_read += read;
		file_pos   += read;
	}

	fclose(file);

	suffix_array->global_byte_end_idx = file_pos;

	clock_gettime(CLOCK_MONOTONIC, &end_time);
	double elapsed = (end_time.tv_sec - start_time.tv_sec) + 
					 (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
	printf("Time to read and parse CSV: %f\n", elapsed);
	printf("Mb/s: %f\n", ((double)bytes_read / (1024 * 1024)) / elapsed);


	suffix_array->n = suffix_array_mapping.buffer_idx;
	suffix_array->suffix_array = (uint32_t*)malloc(suffix_array->n * sizeof(uint32_t));

	construct_truncated_suffix_array(text.buffer, suffix_array);
	free_buffer_c(&text);

	// Map suffix array quoted bits before remapping suffix array.
	buffer_bit suffix_array_is_quoted_copy;
	init_buffer_bit(&suffix_array_is_quoted_copy, suffix_array->n);
	copy_buffer_bit(suffix_array->is_quoted_bitflag, &suffix_array_is_quoted_copy);

	for (uint32_t i = 0; i < suffix_array->n; ++i) {
		uint32_t idx = suffix_array->suffix_array[i];
		set_buffer_bit(
				suffix_array->is_quoted_bitflag, 
				i, 
				get_buffer_bit(&suffix_array_is_quoted_copy, idx)
				);
	}
	free_buffer_bit(&suffix_array_is_quoted_copy);

	// Remap suffix array indices to original file positions.
	#pragma omp parallel for
	for (uint32_t i = 0; i < suffix_array->n; ++i) {
		suffix_array->suffix_array[i] = suffix_array_mapping.buffer[suffix_array->suffix_array[i]];
	}

	free_buffer_u32(&suffix_array_mapping);
}

#define IMM_MODE (uint8_t)( \
					_SIDD_UBYTE_OPS | \
					_SIDD_CMP_EQUAL_EACH  | \
					_SIDD_NEGATIVE_POLARITY | \
					_SIDD_LEAST_SIGNIFICANT \
					)

inline int strncmp_128(const char* str1, const char* str2, uint64_t n) {
    // Load the strings into __m128i registers
    __m128i v1 = _mm_loadu_si128((__m128i*)str1);
    __m128i v2 = _mm_loadu_si128((__m128i*)str2);

    // Setup imm8 to specify unsigned bytes and equality comparison
    // _SIDD_UBYTE_OPS: Treats the comparison data as unsigned bytes.
    // _SIDD_CMP_EQUAL_EACH: Compares each pair of data, looking for equality.
    // _SIDD_BIT_MASK: The result is a bit mask of the comparison results.
	// _SIDD_NEGATIVE_POLARITY: The comparison is negated.

	const int idx = _mm_cmpistri(v1, v2, IMM_MODE);
	return ((uint8_t)str1[idx] - (uint8_t)str2[idx]) * ((uint64_t)idx < n);
}

inline int strncmp_256(const char* str1, const char* str2, uint64_t n) {
    // Load the strings into __m256i registers
	__m256i v1 = _mm256_loadu_si256((__m256i*)str1);
	__m256i v2 = _mm256_loadu_si256((__m256i*)str2);

	// Do a bytewise comp
	__m256i cmp = _mm256_cmpeq_epi8(v1, v2);
	int mask = ~_mm256_movemask_epi8(cmp);

	if (mask == 0) {
		return 0;
	}
	const int idx = __builtin_ctz(mask);

	return ((uint8_t)str1[idx] - (uint8_t)str2[idx]) * ((uint64_t)idx < n);
}

void recursive_bucket_sort(
	const char* str,
	uint32_t* suffix_array,
	uint32_t* temp_suffix_array,
	uint32_t overall_size,
	uint32_t current_size,
	int max_depth,
	int current_depth
) {
	// Base case
	if (current_depth == max_depth) {
		return;
	}

	if (current_size < 32) {
		// Do insertion sort
		for (size_t i = 1; i < current_size; ++i) {
			size_t j = i;
			// while (j > 0 && strncmp(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
			while (j > 0 && strncmp_128(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
				swap(suffix_array[j], suffix_array[j - 1]);
				--j;
			}
		}
		return;
	}

	// Do bucket sort of first char. Then in each bucket, skip
	// current_depth chars and do bucket sort of next char.

	uint32_t _buckets[NUM_BUCKETS] 	     = {0};
	uint32_t _bucket_starts[NUM_BUCKETS] = {0};

	for (uint32_t i = 0; i < current_size; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		++_buckets[char_val];
	}

	uint32_t offset = 0;
	for (uint32_t i = 0; i < NUM_BUCKETS; ++i) {
		_bucket_starts[i] = offset;
		offset += _buckets[i];
	}

	memcpy(temp_suffix_array, suffix_array, current_size * sizeof(uint32_t));

	for (uint32_t i = 0; i < current_size; ++i) {
		uint32_t char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		temp_suffix_array[_bucket_starts[char_val]] = suffix_array[i];
		++_bucket_starts[char_val];
	}

	memcpy(suffix_array, temp_suffix_array, current_size * sizeof(uint32_t));

	// Recalculate bucket starts
	offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		_bucket_starts[i] = offset;
		offset += _buckets[i];
	}

	// Recursively sort each bucket
	if (current_depth == 0) {
		// #pragma omp parallel for schedule(guided)
		#pragma omp parallel for schedule(dynamic)
		for (uint32_t i = 0; i < 255; ++i) {
			if ((char)i == '\n') continue;

			uint32_t bucket_start = _bucket_starts[i];
			uint32_t bucket_end   = _bucket_starts[i + 1];
			if (bucket_end - bucket_start <= 1) {
				continue;
			}

			recursive_bucket_sort(
				str,
				suffix_array + bucket_start,
				temp_suffix_array + bucket_start,
				overall_size,
				bucket_end - bucket_start,
				max_depth,
				current_depth + 1
			);
		}
		return;
	}

	for (uint32_t i = 0; i < 255; ++i) {
		if ((char)i == '\n') continue;

		int bucket_start = _bucket_starts[i];
		int bucket_end   = _bucket_starts[i + 1];
		if (bucket_end - bucket_start <= 1) {
			continue;
		}

		recursive_bucket_sort(
			str,
			suffix_array + bucket_start,
			temp_suffix_array + bucket_start,
			overall_size,
			bucket_end - bucket_start,
			max_depth,
			current_depth + 1
		);
	}
}

void read_text_into_buffer(
	const char* filename,
	char** buffer,
	uint64_t* buffer_size
) {
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	fseek(file, 0, SEEK_END);
	*buffer_size = ftell(file);
	rewind(file);

	printf("File size: %lu MB\n", *buffer_size / (1024 * 1024));

	*buffer = (char*)malloc(*buffer_size);

	fread(*buffer, 1, *buffer_size, file);
	fclose(file);
}

void construct_truncated_suffix_array(const char* str, SuffixArray* suffix_array) {
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	uint32_t _max_suffix_length = min(suffix_array->max_suffix_length, suffix_array->n);
	uint32_t* temp_suffix_array = (uint32_t*)malloc(suffix_array->n * sizeof(uint32_t));

	// Construct simple suffix array from str
	#pragma omp parallel for
	for (uint64_t i = 0; i < suffix_array->n; ++i) {
		suffix_array->suffix_array[i] = i;
	}

	recursive_bucket_sort(
		str,
		suffix_array->suffix_array,
		temp_suffix_array,
		suffix_array->n,
		suffix_array->n,
		_max_suffix_length,
		0
	);

	free(temp_suffix_array);

	clock_gettime(CLOCK_MONOTONIC, &end_time);
	double elapsed = (end_time.tv_sec - start_time.tv_sec) + 
					 (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
	printf("Time to construct suffix array: %f\n", elapsed);
}


pair_u32 get_substring_positions(
    const char* str,
    const SuffixArray* suffix_array,
    const char* substring
) {
    uint32_t m 	   = strlen(substring);
    uint32_t first = 0;
    uint32_t last  = suffix_array->n - 1;

    uint32_t start = UINT32_MAX;
	uint32_t end   = UINT32_MAX;

	uint32_t cmp_length = min(m, suffix_array->max_suffix_length);

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        uint32_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array->suffix_array[mid], substring, cmp_length) < 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, cmp_length) < 0) {
            first = mid + 1;
        }
        else {
            last = mid - 1;
            start = mid;
        }
    }

    if (start == UINT32_MAX) {
		pair_u32 result = {UINT32_MAX, UINT32_MAX};
        return result;
    }

    // Reset for searching the last occurrence
    first = 0;
	last  = suffix_array->n - 1;
    while (first <= last) {
        uint32_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array->suffix_array[mid], substring, cmp_length) > 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, cmp_length) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end = mid;
        }
    }

	pair_u32 result = {start, end};
	return result;
}

pair_u32 get_substring_positions_file(
	FILE* file,
	const SuffixArray* suffix_array,
    const char* substring
) {
    uint32_t m 	   = strlen(substring);
    uint32_t first = 0;
    uint32_t last  = suffix_array->n - 1;

    uint32_t start = UINT32_MAX;
	uint32_t end   = UINT32_MAX;

	uint32_t cmp_length = min(m, suffix_array->max_suffix_length);

	char* line = (char*)malloc(cmp_length * sizeof(char));

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        uint32_t mid = (first + last) / 2;
		fseek(
				file,
				suffix_array->global_byte_start_idx + (uint64_t)suffix_array->suffix_array[mid],
				SEEK_SET
				);
		fread(line, 1, cmp_length, file);

		for (uint32_t i = 0; i < cmp_length; ++i) {
			line[i] = tolower(line[i]);
		}

        if (strncmp(line, substring, cmp_length) < 0) {
        // if (strncmp_128(line, substring, cmp_length) < 0) {
            first = mid + 1;
        }
        else {
            last  = mid - 1;
            start = mid;
        }
    }

    if (start == UINT32_MAX) {
		pair_u32 result = {UINT32_MAX, UINT32_MAX};
        return result;
    }

    // Reset for searching the last occurrence
    first = 0;
	last  = suffix_array->n - 1;
    while (first <= last) {
        uint32_t mid = (first + last) / 2;
		fseek(
				file,
				suffix_array->global_byte_start_idx + (uint64_t)suffix_array->suffix_array[mid],
				SEEK_SET
				);
		fread(line, 1, cmp_length, file);

		for (uint32_t i = 0; i < cmp_length; ++i) {
			line[i] = tolower(line[i]);
		}

        if (strncmp(line, substring, cmp_length) > 0) {
        // if (strncmp_128(line, substring, cmp_length) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end   = mid;
        }
    }

	free(line);
	
	pair_u32 result = {start, end};
	return result;
}


uint32_t get_matching_records(
	const char* str,
	const SuffixArray* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records
) {
	pair_u32 match_idxs = get_substring_positions(
			str, 
			suffix_array, 
			substring
			);

	if (match_idxs.first == UINT32_MAX) {
		return 0;
	}

	uint32_t num_matches = min(k, match_idxs.second - match_idxs.first + 1);

	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		// Go to the original index and iterate backwards until newline.
		uint32_t offset = suffix_array->suffix_array[i];
		while (str[offset] != '\n') --offset;
		++offset;

		uint32_t match_length = 0;
		while (1) {
			char c = str[offset++];
			if (c == '\\') {
				++offset;
				match_length += 2;
				continue;
			}
			if (c == '\n') {
				break;
			}
			++match_length;
		}

		char* record = (char*)malloc(match_length + 1);
		memcpy(record, str + offset, match_length);

		matching_records[i - match_idxs.first] = record;
	}

	return num_matches;
}


static inline uint32_t rfc4180_seek_backward_newline(
	FILE* file,
	uint64_t idx,
	uint8_t in_quotes,
	uint64_t total_bytes_file
) {
	uint64_t num_chars_back = 0;
	char 	 line[DEFAULT_PAGE_SIZE];

	uint8_t init_in_quotes = in_quotes;

	while (1) {
		if (num_chars_back > 16384) {
			printf("Error: Newline not found behind.\n");

			char small_line[128];
			fseek(file, idx - 64, SEEK_SET);
			fread(small_line, 1, 127, file);
			small_line[127] = '\0';
			printf("Line: %s\n", small_line);

			exit(1);
		}

		uint64_t seek_point = 0;
		if (idx <= (num_chars_back + DEFAULT_PAGE_SIZE)) {
			seek_point = 0;
		} else {
			seek_point = idx - (num_chars_back + DEFAULT_PAGE_SIZE);
		}

		assert(seek_point <= idx);

		fseek(file, seek_point, SEEK_SET);
		fread(line, 1, DEFAULT_PAGE_SIZE, file);

		uint32_t start_idx = min(DEFAULT_PAGE_SIZE - 1, idx - seek_point);

		for (uint32_t char_idx = start_idx; char_idx > 0; --char_idx) {
			if (line[char_idx] == '"') {
				if (in_quotes) {
					if (line[char_idx - 1] == '"') {
						--char_idx;
						continue;
					}
					in_quotes = !in_quotes;
				} else {
					in_quotes = !in_quotes;
				}
			}

			if (line[char_idx] == '\n' && !in_quotes) {
				return (uint32_t)(num_chars_back + start_idx - char_idx - 1);
			}
		}

		if (seek_point == 0) break;
		num_chars_back += DEFAULT_PAGE_SIZE;
	}

	return UINT32_MAX;
}

static inline uint32_t rfc4180_seek_forward_newline(
	FILE* file,
	uint64_t idx,
	uint8_t in_quotes,
	uint64_t total_bytes_file
) {
	uint32_t num_chars_ahead = 0;
	char 	 line[DEFAULT_PAGE_SIZE];

	while (1) {
		if (num_chars_ahead > 16777216) {
			printf("Error: Newline not found ahead.\n");
			printf("In quotes: %d\n", in_quotes);
			exit(1);
		}

		uint64_t seek_point = min(total_bytes_file, idx + num_chars_ahead);
		fseek(file, seek_point, SEEK_SET);
		fread(line, 1, DEFAULT_PAGE_SIZE, file);

		for (uint32_t char_idx = 0; char_idx < DEFAULT_PAGE_SIZE - 1; ++char_idx) {
			if (line[char_idx] == '"') {
				if (in_quotes) {
					if (line[char_idx + 1] == '"') {
						++char_idx;
						continue;
					}
					in_quotes = !in_quotes;
				} else {
					in_quotes = !in_quotes;
				}
			}

			if (line[char_idx] == '\n' && !in_quotes) {
				return num_chars_ahead + char_idx - 1;
			}
		}

		if (seek_point == total_bytes_file) break;
		num_chars_ahead += DEFAULT_PAGE_SIZE;
	}

	return UINT32_MAX;
}


uint32_t get_matching_records_file(
	const char* filename,
	const SuffixArray* suffix_array,
	const char* substring,
	uint32_t k,
	char** matching_records
) {
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	pair_u32 match_idxs = get_substring_positions_file(
			file,
			suffix_array,
			substring
			);

	if (match_idxs.first == UINT32_MAX) {
		return 0;
	}

	// uint32_t num_matches = min(k, match_idxs.second - match_idxs.first + 1);
	uint32_t num_matches = min(k, match_idxs.second - match_idxs.first - 1);


	fseek(file, 0, SEEK_END);
	uint64_t total_bytes_file = ftell(file);

	assert(suffix_array->global_byte_start_idx < total_bytes_file);
	assert(suffix_array->global_byte_end_idx  == total_bytes_file);
	assert(suffix_array->is_quoted_bitflag->buffer_bit_idx == suffix_array->n);

	printf("Num matches: %d\n", num_matches);

	uint32_t match_idx = 0;
	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		uint8_t is_quoted 	 = get_buffer_bit(suffix_array->is_quoted_bitflag, i);
		uint64_t byte_offset = (uint64_t)suffix_array->suffix_array[i] + 
										 suffix_array->global_byte_start_idx;

		const uint32_t num = 25;
		for (uint32_t j = i - num; j < i + num; ++j) {
			get_buffer_bit(suffix_array->is_quoted_bitflag, j);
		}

		uint32_t num_chars_back  = rfc4180_seek_backward_newline(
				file, 
				byte_offset, 
				is_quoted,
				total_bytes_file
				);
		uint32_t num_chars_ahead = rfc4180_seek_forward_newline(
				file, 
				byte_offset, 
				is_quoted, 
				total_bytes_file
				);

		char* record = (char*)malloc(num_chars_back + num_chars_ahead + 1);
		fseek(file, byte_offset - num_chars_back, SEEK_SET);
		fread(record, 1, num_chars_back + num_chars_ahead, file);
		record[num_chars_back + num_chars_ahead] = '\0';

		matching_records[match_idx++] = record;
	}
	fclose(file);

	return num_matches;
}
