#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#include <vector>
#include <algorithm>

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
			// while (j > 0 && strncmp(str + suffix_array[j], str + suffix_array[j - 1], max_depth) < 0) {
			while (j > 0 && strncmp(str + suffix_array[j] + current_depth, str + suffix_array[j - 1] + current_depth, max_depth) < 0) {
				std::swap(suffix_array[j], suffix_array[j - 1]);
				--j;
			}
		}
		return;
	}

	// Do bucket sort of first char. Then in each bucket, skip
	// current_depth chars and do bucket sort of next char.
	constexpr int NUM_BUCKETS 	   		= 256;
	uint32_t buckets[NUM_BUCKETS] 	    = {0};
	uint32_t bucket_starts[NUM_BUCKETS] = {0};

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		++buckets[char_val];
	}

	uint64_t offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	// temp
	memcpy(temp_suffix_array, suffix_array, n * sizeof(uint32_t));

	for (uint64_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = (uint8_t)str[char_idx];
		temp_suffix_array[bucket_starts[char_val]] = suffix_array[i];
		++bucket_starts[char_val];
	}

	memcpy(suffix_array, temp_suffix_array, n * sizeof(uint32_t));

	// Recalculate bucket starts
	offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	// Recursively sort each bucket
	if (current_depth == 0) {
		#pragma omp parallel for schedule(static, 1)
		for (int i = 0; i < NUM_BUCKETS - 1; ++i) {
			int bucket_start = bucket_starts[i];
			int bucket_end   = bucket_starts[i + 1];
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
	for (int i = 0; i < NUM_BUCKETS - 1; ++i) {
		int bucket_start = bucket_starts[i];
		int bucket_end   = bucket_starts[i + 1];
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
	uint32_t max_suffix_length
) {
	max_suffix_length = std::min(max_suffix_length, (uint32_t)n);
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

	// Construct simple suffix array from str
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

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

}


std::vector<uint32_t> get_substring_positions(
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
            first = mid + 1;
        }
        else {
            last = mid - 1;
            start = mid;
        }
    }

    if (start == -1) {
        return std::vector<uint32_t>();
    }

    // Reset for searching the last occurrence
    first = 0, last = n - 1;
    while (first <= last) {
        int64_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array[mid], substring, m) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end = mid;
        }
    }

	std::vector<uint32_t> matches;
	matches.reserve(end - start + 1);
	for (int i = start; i <= end; ++i) {
		matches.push_back(suffix_array[i]);
	}
	return matches;
}

std::vector<uint32_t> get_matching_indices(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	const char* substring,
	const uint32_t* row_offsets,
	uint32_t num_rows,
	int k
) {
	// First get matches indices from the suffix array
	// then convert these to row indices again with binary search.
	std::vector<uint32_t> matches = get_substring_positions(str, suffix_array, n, substring);
	printf("Num matches: %lu\n", matches.size());

	size_t num_matches = std::min((size_t)k, matches.size());

	std::vector<uint32_t> row_indices;
	row_indices.reserve(num_matches);

	for (uint32_t i = 0; i < num_matches; ++i) {
		uint32_t match = matches[i];
		uint32_t row = 0;
		uint32_t first = 0;
		uint32_t last = num_rows - 1;
		while (first <= last) {
			uint32_t mid = (first + last) / 2;
			if (row_offsets[mid] <= match) {
				row = mid;
				first = mid + 1;
			}
			else {
				last = mid - 1;
			}
		}

		row_indices.push_back(row);
	}

	return row_indices;
}
