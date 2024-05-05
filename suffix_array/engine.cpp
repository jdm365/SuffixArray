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


void construct_truncated_suffix_array_from_csv(
	const char* csv_file,
	uint32_t column_idx,
	uint32_t* suffix_array,
	uint32_t max_suffix_length
) {
	auto start = std::chrono::high_resolution_clock::now();

	// Read and parse the CSV file.
	FILE* file = fopen(csv_file, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}

	char*    line = NULL;
	size_t   len = 0;
	ssize_t  read;
	uint64_t file_pos = 0;
	uint64_t suffix_array_size = 0;

	uint64_t line_0_size = getline(&line, &len, file);
	rewind(file);

	fseek(file, 0, SEEK_END);
	uint64_t num_lines = ftell(file) / line_0_size;
	rewind(file);

	std::vector<char> text;
	text.reserve(num_lines * max_suffix_length);

	suffix_array = (uint32_t*)malloc(num_lines * max_suffix_length * sizeof(uint32_t));

	while ((read = getline(&line, &len, file)) != -1) {
		uint32_t char_idx = 0;
		uint32_t col_idx  = 0;

		while (col_idx < column_idx) {
			if (line[char_idx] == ',') {
				++col_idx;
			}
			++char_idx;
		}

		while (line[char_idx] != '\n' && line[char_idx] != '\0' && line[char_idx] != ',') {
			text.push_back(line[char_idx]);
			suffix_array[suffix_array_size++] = file_pos + char_idx;
			++char_idx;
		}

		text.push_back('\n');
		suffix_array[suffix_array_size++] = file_pos + char_idx;

		file_pos += read;
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	printf("Time to read and parse CSV: %f\n", elapsed.count());

	suffix_array = (uint32_t*)realloc(suffix_array, suffix_array_size * sizeof(uint32_t));

	_construct_truncated_suffix_array_preset(
		text.data(),
		suffix_array,
		suffix_array_size,
		max_suffix_length
	);

	fclose(file);
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
	int current_depth
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
			// while (j > 0 && strncmp_256(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
				std::swap(suffix_array[j], suffix_array[j - 1]);
				--j;
			}
		}
		return;
	}

	// Do bucket sort of first char. Then in each bucket, skip
	// current_depth chars and do bucket sort of next char.

	/*
	const int thread_id = omp_get_thread_num();
	const int thread_offset = thread_id * 256;
	for (size_t i = 0; i < 256; ++i) {
		buckets[thread_offset + i] = 0;
	}

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		++buckets[thread_offset + char_val];
	}

	uint64_t offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		bucket_starts[i + thread_offset] = offset;
		offset += buckets[i + thread_offset];
	}

	// temp
	memcpy(temp_suffix_array, suffix_array, n * sizeof(uint32_t));

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		temp_suffix_array[bucket_starts[char_val + thread_offset]] = suffix_array[i];
		++bucket_starts[char_val + thread_offset];
	}

	memcpy(suffix_array, temp_suffix_array, n * sizeof(uint32_t));
	*/
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
				current_depth + 1
			);
		}
		return;
	}

	for (int i = 0; i < 255; ++i) {
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

	fflush(stdout);
}

void construct_truncated_suffix_array(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	uint32_t max_suffix_length,
	bool use_index_array
) {
	max_suffix_length = std::min(max_suffix_length, (uint32_t)n);
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

	// Construct simple suffix array from str
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

	const int num_threads = omp_get_max_threads();

	// alignas(64) uint32_t* buckets 		= (uint32_t*)malloc(num_threads * 256 * sizeof(uint32_t));
	// alignas(64) uint32_t* bucket_starts = (uint32_t*)malloc(num_threads * 256 * sizeof(uint32_t));

	recursive_bucket_sort(
		str,
		suffix_array,
		temp_suffix_array,
		n,
		n,
		max_suffix_length,
		0
	);

	free(temp_suffix_array);

	// free(buckets);
	// free(bucket_starts);
}

void _construct_truncated_suffix_array_preset(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	uint32_t max_suffix_length
) {
	printf("Constructing suffix array of size: %lu\n", n);

	max_suffix_length = std::min(max_suffix_length, (uint32_t)n);
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));
	memcpy(temp_suffix_array, suffix_array, n * sizeof(uint32_t));

	// const int num_threads = omp_get_max_threads();
	// alignas(64) uint32_t* buckets 		= (uint32_t*)malloc(num_threads * 256 * sizeof(uint32_t));
	// alignas(64) uint32_t* bucket_starts = (uint32_t*)malloc(num_threads * 256 * sizeof(uint32_t));

	recursive_bucket_sort(
		str,
		suffix_array,
		temp_suffix_array,
		// buckets,
		// bucket_starts,
		n,
		n,
		max_suffix_length,
		0
	);

	free(temp_suffix_array);

	// free(buckets);
	// free(bucket_starts);
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

std::vector<uint32_t> get_matching_indices(
	const char* str,
	uint32_t* suffix_array,
	uint32_t* suffix_array_idxs,
	uint64_t n,
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
		return std::vector<uint32_t>();
	}

	size_t num_matches = std::min((size_t)k, (size_t)(match_idxs.second - match_idxs.first + 1));

	std::vector<uint32_t> matches;
	matches.reserve(num_matches);

	robin_hood::unordered_flat_set<uint32_t> match_set;

	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		if (match_set.insert(suffix_array_idxs[suffix_array[i]]).second) {
			matches.push_back(suffix_array_idxs[suffix_array[i]]);
		}
	}
	return matches;
}


std::vector<uint32_t> get_matching_indices_no_idxs(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
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
		return std::vector<uint32_t>();
	}

	size_t num_matches = std::min((size_t)k, (size_t)(match_idxs.second - match_idxs.first + 1));

	std::vector<uint32_t> matches;
	matches.reserve(num_matches);

	robin_hood::unordered_flat_set<uint32_t> match_set;

	for (uint32_t i = match_idxs.first; i < match_idxs.first + num_matches; ++i) {
		// Go to the original index and iterate backwards until newline.
		uint32_t offset = suffix_array[i];
		while (str[offset] != '\n') --offset;
		++offset;

		if (match_set.insert(offset).second) {
			matches.push_back(offset);
		}
	}
	return matches;
}
