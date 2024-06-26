#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>

#include <immintrin.h>
#include "robin_hood.h"

#include "engine.h"


inline uint32_t min(uint32_t a, uint32_t b) {
	return a < b ? a : b;
}

inline uint32_t max(uint32_t a, uint32_t b) {
	return a > b ? a : b;
}

inline uint64_t min(uint64_t a, uint64_t b) {
	return a < b ? a : b;
}

inline uint64_t max(uint64_t a, uint64_t b) {
	return a > b ? a : b;
}

inline void parse_line_v2(
		const char* line,
		std::vector<char>& text,
		std::vector<uint32_t>& suffix_array_mapping,
		uint64_t file_pos,
		uint32_t column_idx
		) {

	uint32_t char_idx = 0;
	uint32_t col_idx  = 0;

	uint32_t line_length = strlen(line);
	bool first_in_col = true;

	__m128i comma_mask   = _mm_set1_epi8(',');
	__m128i newline_mask = _mm_set1_epi8('\n');
	__m128i quote_mask   = _mm_set1_epi8('"');
	__m128i escape_mask = _mm_set1_epi8('\\');

	__m128i zero = _mm_setzero_si128();

	while (col_idx < column_idx) {
		__m128i line_v = _mm_loadu_si128((__m128i*)(line + char_idx));
		__m128i comma_cmp = _mm_cmpeq_epi8(line_v, comma_mask);
		__m128i newline_cmp = _mm_cmpeq_epi8(line_v, newline_mask);
		__m128i quote_cmp = _mm_cmpeq_epi8(line_v, quote_mask);
		__m128i escape_cmp = _mm_cmpeq_epi8(line_v, escape_mask);

		__m128i comma_masked = _mm_and_si128(comma_cmp, _mm_cmplt_epi8(zero, line_v));
		__m128i newline_masked = _mm_and_si128(newline_cmp, _mm_cmplt_epi8(zero, line_v));
		__m128i quote_masked = _mm_and_si128(quote_cmp, _mm_cmplt_epi8(zero, line_v));
		__m128i escape_masked = _mm_and_si128(escape_cmp, _mm_cmplt_epi8(zero, line_v));

		__m128i mask = _mm_or_si128(
			_mm_or_si128(
				_mm_or_si128(comma_masked, newline_masked),
				quote_masked
			),
			escape_masked
		);

		uint32_t mask_int = _mm_movemask_epi8(mask);
		if (mask_int != 0) {
			uint32_t idx = __builtin_ctz(mask_int);
			char_idx += idx + 1;
			++col_idx;
			first_in_col = true;
			continue;
		}

		char_idx += 16;
	}

	first_in_col = true;

	for (; char_idx < line_length;) {
		__m128i line_v      = _mm_loadu_si128((__m128i*)(line + char_idx));
		__m128i comma_cmp   = _mm_cmpeq_epi8(line_v, comma_mask);
		__m128i newline_cmp = _mm_cmpeq_epi8(line_v, newline_mask);
		__m128i quote_cmp   = _mm_cmpeq_epi8(line_v, quote_mask);
		__m128i escape_cmp  = _mm_cmpeq_epi8(line_v, escape_mask);

		__m128i comma_masked   = _mm_and_si128(comma_cmp,   _mm_cmplt_epi8(zero, line_v));
		__m128i newline_masked = _mm_and_si128(newline_cmp, _mm_cmplt_epi8(zero, line_v));
		__m128i quote_masked   = _mm_and_si128(quote_cmp,   _mm_cmplt_epi8(zero, line_v));
		__m128i escape_masked  = _mm_and_si128(escape_cmp,  _mm_cmplt_epi8(zero, line_v));

		__m128i mask = _mm_or_si128(
			_mm_or_si128(
				_mm_or_si128(comma_masked, newline_masked),
				quote_masked
			),
			escape_masked
		);

		uint32_t mask_int = _mm_movemask_epi8(mask);
		if (mask_int != 0) {
			uint32_t idx = __builtin_ctz(mask_int);
			char_idx += idx + 1;
			break;
		}

		__m128i text_v = _mm_and_si128(line_v, _mm_set1_epi8(0xDF));
		__m128i suffix_v = _mm_set1_epi32(file_pos + char_idx);

		_mm_storeu_si128((__m128i*)(text.data() + text.size()), text_v);
		_mm_storeu_si128((__m128i*)(suffix_array_mapping.data() + suffix_array_mapping.size()), suffix_v);

		text.push_back(tolower(line[char_idx]));
		suffix_array_mapping.push_back(file_pos + char_idx);
		++char_idx;

		first_in_col = false;
	}

	/*
	while (line[char_idx] != '\n' && line[char_idx] != '\0' && line[char_idx] != ',') {
		if (line[char_idx] == '\\') {
			++char_idx;
			text.push_back(tolower(line[char_idx]));
			suffix_array_mapping.push_back(file_pos + char_idx);
			++char_idx;

			first_in_col = false;
			continue;
		}
		if (line[char_idx] == '"' && first_in_col) {
			first_in_col = false;

			++char_idx;
			while (line[char_idx] != '"') {
				if (char_idx >= line_length) {
					return;
				}
				text.push_back(tolower(line[char_idx]));
				suffix_array_mapping.push_back(file_pos + char_idx);
				++char_idx;
			}
			++char_idx;
			continue;
		}

		text.push_back(tolower(line[char_idx]));
		suffix_array_mapping.push_back(file_pos + char_idx);
		++char_idx;

		first_in_col = false;
	}
	*/

	text.push_back('\n');
	suffix_array_mapping.push_back(file_pos + char_idx);
}

inline void parse_line(
		const char* line,
		std::vector<char>& text,
		std::vector<uint32_t>& suffix_array_mapping,
		uint64_t file_pos,
		uint32_t column_idx
		) {

	uint32_t char_idx = 0;
	uint32_t col_idx  = 0;

	uint32_t line_length = strlen(line);
	bool first_in_col = true;

	while (col_idx < column_idx) {
		if (line[char_idx] == '\\') {
			char_idx += 2;
			first_in_col = false;
			continue;
		}

		if (line[char_idx] == '"' && first_in_col) {
			first_in_col = false;
			++char_idx;
			while (line[char_idx] != '"') {
				++char_idx;
				if (char_idx >= line_length) {
					/*
					printf("Error: Invalid CSV format\n");
					printf("Line: %s\n\n", line);
					printf("File position: %lu\n", file_pos);
					*/
					return;
					// exit(1);
				}
			}
			++char_idx;
		}
		if (line[char_idx] == ',') {
			++col_idx;
			first_in_col = true;
		} else {
			first_in_col = false;
		}
		++char_idx;
	}

	first_in_col = true;

	while (line[char_idx] != '\n' && line[char_idx] != '\0' && line[char_idx] != ',') {
		if (line[char_idx] == '\\') {
			++char_idx;
			text.push_back(tolower(line[char_idx]));
			suffix_array_mapping.push_back(file_pos + char_idx);
			++char_idx;

			first_in_col = false;
			continue;
		}
		if (line[char_idx] == '"' && first_in_col) {
			first_in_col = false;

			++char_idx;
			while (line[char_idx] != '"') {
				if (char_idx >= line_length) {
					/*
					printf("Error: Invalid CSV format\n");
					printf("Line: %s\n\n", line);
					printf("File position: %lu\n", file_pos);
					*/
					return;
					// exit(1);
				}
				text.push_back(tolower(line[char_idx]));
				suffix_array_mapping.push_back(file_pos + char_idx);
				++char_idx;
			}
			++char_idx;
			continue;
		}

		text.push_back(tolower(line[char_idx]));
		suffix_array_mapping.push_back(file_pos + char_idx);
		++char_idx;

		first_in_col = false;
	}

	text.push_back('\n');
	suffix_array_mapping.push_back(file_pos + char_idx);
}


void construct_truncated_suffix_array_from_csv_partitioned(
	const char* csv_file,
	uint32_t column_idx,
	std::vector<uint32_t>& suffix_array,
	uint32_t* suffix_array_size,
	uint32_t max_suffix_length,
	uint64_t start_idx,
	uint64_t& end_idx,
	bool skip_newline
) {
	const uint32_t TWO_GB = (uint32_t)2 * (uint32_t)1024 * (uint32_t)1024 * (uint32_t)1024;

	auto start = std::chrono::high_resolution_clock::now();

	// Read and parse the CSV file.
	FILE* file = fopen(csv_file, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	uint32_t max_line_length = 1024;

	char*    line = NULL;
	uint64_t len = 0;
	ssize_t  read;
	uint64_t file_pos = start_idx;
	uint64_t bytes_read = 0;

	uint64_t line_0_size = getline(&line, &len, file);

	uint64_t num_lines = TWO_GB / line_0_size;
	fseek(file, start_idx, SEEK_SET);

	std::vector<char> text;
	text.reserve(num_lines * max_suffix_length);

	std::vector<uint32_t> suffix_array_mapping;
	suffix_array_mapping.reserve(num_lines * max_suffix_length);

	while (
			((read = getline(&line, &len, file)) != -1)
				&& 
			(bytes_read < TWO_GB)
			) {
		if (read > max_line_length) max_line_length = read;

		parse_line(
			line,
			text,
			suffix_array_mapping,
			// file_pos,
			bytes_read,
			column_idx
		);

		bytes_read += read;
		file_pos   += read;
	}

	fclose(file);

	end_idx = file_pos;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	printf("Time to read and parse CSV: %f\n", elapsed.count());
	printf("Mb/s: %f\n", (double)bytes_read / (1024 * 1024) / elapsed.count());

	*suffix_array_size = text.size();
	suffix_array.resize(*suffix_array_size);

	construct_truncated_suffix_array(
		text.data(),
		suffix_array,
		*suffix_array_size,
		max_suffix_length,
		skip_newline
	);

	// Remap suffix array indices to original file positions.
	#pragma omp parallel for
	for (uint32_t i = 0; i < *suffix_array_size; ++i) {
		suffix_array[i] = suffix_array_mapping[suffix_array[i]];
	}
}

static inline int strncmp_128(const char* str1, const char* str2, uint64_t n) {
    // Load the strings into __m128i registers
    __m128i v1 = _mm_loadu_si128((__m128i*)str1);
    __m128i v2 = _mm_loadu_si128((__m128i*)str2);

    // Setup imm8 to specify unsigned bytes and equality comparison
    // _SIDD_UBYTE_OPS: Treats the comparison data as unsigned bytes.
    // _SIDD_CMP_EQUAL_EACH: Compares each pair of data, looking for equality.
    // _SIDD_BIT_MASK: The result is a bit mask of the comparison results.
	// _SIDD_NEGATIVE_POLARITY: The comparison is negated.
	constexpr const uint8_t mode =
            _SIDD_UBYTE_OPS |
            _SIDD_CMP_EQUAL_EACH |
            _SIDD_NEGATIVE_POLARITY |
            _SIDD_LEAST_SIGNIFICANT;


	const int idx = _mm_cmpistri(v1, v2, mode);
	return ((uint8_t)str1[idx] - (uint8_t)str2[idx]) * ((uint64_t)idx < n);
}

static inline int strncmp_256(const char* str1, const char* str2, uint64_t n) {
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
	int string_length,
	uint64_t n,
	int max_depth,
	int current_depth,
	bool skip_newline
) {
	// Base case
	if (current_depth == max_depth) {
		return;
	}

	if (n < 32) {
		// Do insertion sort
		for (uint64_t i = 1; i < n; ++i) {
			int j = i;
			// while (j > 0 && strncmp(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
			while (j > 0 && strncmp_128(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
				std::swap(suffix_array[j], suffix_array[j - 1]);
				--j;
			}
		}
		return;
	}

	// Do bucket sort of first char. Then in each bucket, skip
	// current_depth chars and do bucket sort of next char.

	constexpr int NUM_BUCKETS 	   		= 256;
	uint32_t _buckets[NUM_BUCKETS] 	    = {0};
	uint32_t _bucket_starts[NUM_BUCKETS] = {0};

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		++_buckets[char_val];
	}

	uint64_t offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		_bucket_starts[i] = offset;
		offset += _buckets[i];
	}

	// temp
	memcpy(temp_suffix_array, suffix_array, n * sizeof(uint32_t));

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		temp_suffix_array[_bucket_starts[char_val]] = suffix_array[i];
		++_bucket_starts[char_val];
	}

	memcpy(suffix_array, temp_suffix_array, n * sizeof(uint32_t));

	// Recalculate bucket starts
	offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		_bucket_starts[i] = offset;
		offset += _buckets[i];
	}

	// Recursively sort each bucket
	if (current_depth == 0) {
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < 255; ++i) {

			// Don't sort newline characters
			if (skip_newline && (char)i == '\n') {
				continue;
			}

			int bucket_start = _bucket_starts[i];
			int bucket_end   = _bucket_starts[i + 1];
			if (bucket_end - bucket_start <= 1) {
				continue;
			}

			recursive_bucket_sort(
				str,
				suffix_array + bucket_start,
				temp_suffix_array + bucket_start,
				string_length,
				bucket_end - bucket_start,
				max_depth,
				current_depth + 1,
				skip_newline
			);
		}
		return;
	}

	for (int i = 0; i < 255; ++i) {
		if (skip_newline && (char)i == '\n') {
			continue;
		}

		int bucket_start = _bucket_starts[i];
		int bucket_end   = _bucket_starts[i + 1];
		if (bucket_end - bucket_start <= 1) {
			continue;
		}

		recursive_bucket_sort(
			str,
			suffix_array + bucket_start,
			temp_suffix_array + bucket_start,
			string_length,
			bucket_end - bucket_start,
			max_depth,
			current_depth + 1,
			skip_newline
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

void construct_truncated_suffix_array(
	const char* str,
	std::vector<uint32_t>& suffix_array,
	uint32_t n,
	uint32_t max_suffix_length,
	bool skip_newline
) {

	max_suffix_length = std::min(max_suffix_length, (uint32_t)n);
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

	// Construct simple suffix array from str
	#pragma omp parallel for
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

	recursive_bucket_sort(
		str,
		suffix_array.data(),
		temp_suffix_array,
		n,
		n,
		max_suffix_length,
		0,
		skip_newline
	);

	free(temp_suffix_array);
}


std::pair<uint32_t, uint32_t> get_substring_positions(
    const char* str,
    uint32_t* suffix_array,
    uint64_t n,
    const char* substring
) {
    int64_t m = strlen(substring);
    int64_t first = 0;
    int64_t last = n - 1;
    int64_t start = -1;
	int64_t end = -1;

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        int64_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array[mid], substring, m) < 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, m) < 0) {
            first = mid + 1;
        }
        else {
            last = mid - 1;
            start = mid;
        }
    }

    if (start == -1) {
        return std::make_pair(-1, -1);
    }

    // Reset for searching the last occurrence
    first = 0, last = n - 1;
    while (first <= last) {
        int64_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array[mid], substring, m) > 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, m) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end = mid;
        }
    }

	return std::make_pair(start, end);
}

std::pair<uint64_t, uint64_t> get_substring_positions_file(
	FILE* file,
    uint64_t byte_offset,
    uint32_t* suffix_array,
    uint32_t n,
    const char* substring
) {
    int64_t m = strlen(substring);
    int64_t first = 0;
    int64_t last = n - 1;

    int64_t start = -1;
	int64_t end = -1;

	char* line = NULL;
	line = (char*)malloc(32 * sizeof(char));

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        int64_t mid = (first + last) / 2;
		fseek(file, byte_offset + (uint64_t)suffix_array[mid], SEEK_SET);
		fread(line, 1, 32, file);

		for (int i = 0; i < 32; ++i) {
			line[i] = tolower(line[i]);
		}

        if (strncmp(line, substring, m) < 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, m) < 0) {
            first = mid + 1;
        }
        else {
            last = mid - 1;
            start = mid;
        }
    }

    if (start == -1) {
        return std::make_pair(UINT32_MAX, UINT32_MAX);
    }

    // Reset for searching the last occurrence
    first = 0, last = n - 1;
    while (first <= last) {
        int64_t mid = (first + last) / 2;
		fseek(file, byte_offset + (uint64_t)suffix_array[mid], SEEK_SET);
		fread(line, 1, 32, file);

		for (int i = 0; i < 32; ++i) {
			line[i] = tolower(line[i]);
		}

        if (strncmp(line, substring, m) > 0) {
        // if (strncmp_128(str + suffix_array[mid], substring, m) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end = mid;
        }
    }

	return std::make_pair(start, end);
}


std::vector<std::string> get_matching_records(
	const char* str,
	uint32_t* suffix_array,
	uint32_t n,
	const char* substring,
	int k 
) {
	std::pair<uint32_t, uint32_t> match_idxs = get_substring_positions(
			str, 
			suffix_array, 
			n, 
			substring
			);

	if ((int)match_idxs.first == -1) {
		return std::vector<std::string>();
	}

	uint32_t num_matches = std::min((size_t)k, (size_t)(match_idxs.second - match_idxs.first + 1));

	std::string record;

	std::vector<std::string> matches;
	matches.reserve(num_matches);

	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		// Go to the original index and iterate backwards until newline.
		uint32_t offset = suffix_array[i];
		while (str[offset] != '\n') --offset;
		++offset;

		while (true) {
			char c = str[offset++];
			if (c == '\\') {
				++offset;
				record.push_back(str[offset++]);
				continue;
			}
			if (c == '\n') {
				break;
			}
			record.push_back(c);
		}

		matches.push_back(record);
		record.clear();
	}
	return matches;
}

std::vector<std::string> get_matching_records_file(
	const char* filename,
	uint32_t* suffix_array,
	uint32_t n,
	uint64_t byte_offset,
	const char* substring,
	int k 
) {
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	std::pair<uint64_t, uint64_t> match_idxs = get_substring_positions_file(
			file,
			byte_offset,
			suffix_array, 
			n, 
			substring
			);

	if (match_idxs.first == UINT32_MAX) {
		return std::vector<std::string>();
	}

	size_t num_matches = std::min((size_t)k, (size_t)(match_idxs.second - match_idxs.first + 1));

	std::vector<std::string> records;
	records.reserve(num_matches);

	robin_hood::unordered_flat_set<uint32_t> match_set;

	char line[1024];

	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		// Go to the original index and iterate backwards until newline.
		uint64_t offset = suffix_array[i] + byte_offset;
		uint32_t newline_pos = 0;

		fseek(file, max(offset - 512, 0), SEEK_SET);
		fread(line, 1, 1024, file);

		for (uint32_t j = 0; j < 512; ++j) {
			if (line[j] == '\\') {
				j += 2;
				continue;
			}
			if (line[j] == '\n') {
				newline_pos = j;
			}
		}
		offset = offset - 512 + newline_pos + 1;

		fseek(file, offset, SEEK_SET);
		fread(line, 1, 1024, file);

		std::string record;
		uint32_t char_idx = 0;
		while (true) {
			char c = line[char_idx++];
			if (c == '\n') {
				break;
			}
			record.push_back(c);
		}

		records.push_back(record);
	}
	fclose(file);

	return records;
}
