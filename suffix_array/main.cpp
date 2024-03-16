#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

// #include "libsais/include/libsais64.h"
#include "divsufsort.h"

static void print_suffix_array(uint32_t* suffix_array, uint64_t n) {
	for (uint32_t i = 0; i < n; ++i) {
		printf("%d ", suffix_array[i]);
	}
	printf("\n");
}

static void print_suffix_array_chars(const char* str, int* suffix_array, uint64_t n) {
	for (uint32_t i = 0; i < n; ++i) {
		printf("%u ", (uint8_t)str[suffix_array[i]]);
	}
	printf("\n");
}

void two_char_bucket_sort(
		const char* str,
		uint32_t* bucket_starts,
		uint32_t* suffix_array,
		uint64_t n
) {
	constexpr int TWO_CHAR_BUCKETS = 65536;

	uint32_t buckets[TWO_CHAR_BUCKETS] = {0};

	// Do initial bucket sort up to 2 char length (bucket sort)
	for (uint32_t i = 0; i < (int64_t)n - 1; ++i) {
		uint16_t two_char = (uint16_t)str[i] * 256 + (uint16_t)str[i + 1];
		++buckets[two_char];
	}

	uint32_t offset = 0;

	for (size_t i = 0; i < TWO_CHAR_BUCKETS; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	for (int64_t i = 0; i < (int64_t)n - 1; ++i) {
		uint16_t two_char = (uint16_t)str[i] * 256 + (uint16_t)str[i + 1];
		uint64_t new_idx = bucket_starts[two_char];

		suffix_array[new_idx] = i;
		++bucket_starts[two_char];
	}
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

	int offset = 0;
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
		// printf("Bucket start: %d, end: %d\n", bucket_start, bucket_end);

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



void construct_suffix_array(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n
) {
	// Construct simple suffix array from str
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

	// Sort suffix array
	std::sort(
		suffix_array,
		suffix_array + n,
		[&](int a, int b) {
			return strncmp(str + a, str + b, 12) < 0;
		}
	);
}

void construct_truncated_suffix_array(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	int max_suffix_length = -1
) {
	if (max_suffix_length == -1) {
		max_suffix_length = n;
	}
	else {
		max_suffix_length = std::min(max_suffix_length, (int)n);
	}

	// Construct simple suffix array from str
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

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
	return;

	constexpr uint32_t TWO_CHAR_BUCKETS = 256 * 256;

	uint32_t bucket_starts[TWO_CHAR_BUCKETS] = {0};
	two_char_bucket_sort(str, bucket_starts, suffix_array, n);

	uint32_t max_bucket_size = 0;
	for (uint32_t i = 0; i < TWO_CHAR_BUCKETS - 1; ++i) {
		max_bucket_size = std::max(max_bucket_size, bucket_starts[i + 1] - bucket_starts[i]);
	}


	int num_threads = omp_get_max_threads();
	int CHUNKSIZE = 1;


	// Recalculate bucket starts
	uint32_t buckets[TWO_CHAR_BUCKETS] = {0};
	memset(bucket_starts, 0, TWO_CHAR_BUCKETS * sizeof(uint32_t));
	for (uint32_t i = 0; i < (int64_t)n - 1; ++i) {
		uint16_t two_char = (uint16_t)str[i] * 256 + (uint16_t)str[i + 1];
		++buckets[two_char];
	}

	uint32_t offset = 0;

	for (size_t i = 0; i < TWO_CHAR_BUCKETS; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	// #pragma omp parallel for schedule(static, CHUNKSIZE)
	for (uint32_t i = 0; i < TWO_CHAR_BUCKETS - 1; ++i) {
		if (bucket_starts[i + 1] - bucket_starts[i] <= 1) {
			continue;
		}
		// auto start = std::chrono::high_resolution_clock::now();

		recursive_bucket_sort(
				str, 
				suffix_array + bucket_starts[i], 
				temp_suffix_array + bucket_starts[i],
				n, 
				bucket_starts[i + 1] - bucket_starts[i],
				max_suffix_length, 
				2
				);

		/*
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		int microseconds = elapsed.count() * 1000000;

		printf("Bucket %d: %d  %dus\n", i, bucket_starts[i + 1] - bucket_starts[i], microseconds);
		*/
	}

	free(temp_suffix_array);
}

/*

// Forward declarations
void construct_summary_suffix_array(
		int* summary_str,
		int* summary_suffix_array,
		int n,
		int summary_alphabet_size
);
void construct_sais_suffix_array(
	const char* str,
	int* suffix_array,
	int n
);

void induce_sort_l_chars(
	const char* str,
	int* suffix_array,
	uint8_t* types,
	int* buckets,
	int n
) {
	// LMS chars are now placed.
	// Now we can induce sort the L-type chars.
	// We will iterate backwards through the suffix array. 
	// When we encounter an LMS char (not -1), we will check 
	// if the previous char is L-type. If it is, we will place
	// it in the suffix array at the start of its bucket and
	// increment the bucket start.
	int* bucket_starts = (int*)malloc(256 * sizeof(int));
	memset(bucket_starts, 0, 256 * sizeof(int));

	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	for (int i = 0; i < n; ++i) {
		int suffix_idx = suffix_array[i];
		if (suffix_idx < 0) {
			continue;
		}

		uint8_t prev_char_type = types[suffix_idx - 1];
		if (prev_char_type != 1) {
			continue;
		}

		// We found an L-type char.
		int bucket_char = (int)str[suffix_idx - 1];
		int new_idx = bucket_starts[bucket_char];

		suffix_array[new_idx] = suffix_idx - 1;
		++bucket_starts[bucket_char];
	}

	free(bucket_starts);
}

void induce_sort_s_chars(
		const char* str,
		int* suffix_array,
		uint8_t* types,
		int* buckets,
		int n
) {
	// L-type chars are now placed.
	// Now we can induce sort the S-type chars.
	// Using the post-LMS-char placement bucket ends (current non-filled ends)
	// we can again iterate backwards through the suffix array.
	// When we encounter an LMS char (not -1), we will check if the previous
	// char is S-type. If it is, we will place it in the suffix array at the end
	// of its bucket and decrement the bucket end.
	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	memset(bucket_ends, 0, 256 * sizeof(int));

	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	for (int i = n - 1; i >= 0; --i) {
		int suffix_idx = suffix_array[i];
		if (suffix_idx < 1) {
			continue;
		}

		uint8_t prev_char_type = types[suffix_idx - 1];
		if (prev_char_type != 0) {
			continue;
		}

		// We found an S-type char.
		int bucket_char = (int)str[suffix_idx - 1];
		int new_idx = bucket_ends[bucket_char];

		suffix_array[new_idx] = suffix_idx - 1;
		--bucket_ends[bucket_char];
	}

	free(bucket_ends);
}


void induce_sort_l_chars(
	int* str,
	int* suffix_array,
	uint8_t* types,
	int* buckets,
	int n
) {
	// LMS chars are now placed.
	// Now we can induce sort the L-type chars.
	// We will iterate backwards through the suffix array. 
	// When we encounter an LMS char (not -1), we will check 
	// if the previous char is L-type. If it is, we will place
	// it in the suffix array at the start of its bucket and
	// increment the bucket start.
	int* bucket_starts = (int*)malloc(256 * sizeof(int));
	memset(bucket_starts, 0, 256 * sizeof(int));

	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}

	for (int i = 0; i < n; ++i) {
		int suffix_idx = suffix_array[i];
		if (suffix_idx < 0) {
			continue;
		}

		uint8_t prev_char_type = types[suffix_idx - 1];
		if (prev_char_type != 1) {
			continue;
		}

		// We found an L-type char.
		int bucket_char = str[suffix_idx - 1];
		int new_idx = bucket_starts[bucket_char];

		suffix_array[new_idx] = suffix_idx - 1;
		++bucket_starts[bucket_char];
	}

	free(bucket_starts);
}

void induce_sort_s_chars(
		int* str,
		int* suffix_array,
		uint8_t* types,
		int* buckets,
		int n
) {
	// L-type chars are now placed.
	// Now we can induce sort the S-type chars.
	// Using the post-LMS-char placement bucket ends (current non-filled ends)
	// we can again iterate backwards through the suffix array.
	// When we encounter an LMS char (not -1), we will check if the previous
	// char is S-type. If it is, we will place it in the suffix array at the end
	// of its bucket and decrement the bucket end.
	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	memset(bucket_ends, 0, 256 * sizeof(int));

	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	for (int i = n - 1; i >= 0; --i) {
		int suffix_idx = suffix_array[i];
		if (suffix_idx < 1) {
			continue;
		}

		uint8_t prev_char_type = types[suffix_idx - 1];
		if (prev_char_type != 0) {
			continue;
		}

		// We found an S-type char.
		int bucket_char = str[suffix_idx - 1];
		int new_idx = bucket_ends[bucket_char];

		suffix_array[new_idx] = suffix_idx - 1;
		--bucket_ends[bucket_char];
	}

	free(bucket_ends);
}

void accurate_lms_sort(
		const char* str,
		int* buckets,
		uint8_t* types,
		int* summary_suffix_array,
		int* summary_offsets,
		int* final_suffix_array,
		int n
) {
	memset(final_suffix_array, -1, n * sizeof(int));

	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	int  offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	for (int i = n - 1; i >= 0; --i) {
		int str_idx = summary_offsets[summary_suffix_array[i]];

		int char_idx = (int)str[str_idx];
		final_suffix_array[bucket_ends[char_idx]] = str_idx;
		--bucket_ends[char_idx];
	}

	final_suffix_array[0] = n - 1;
}

void accurate_lms_sort(
		int* str,
		int* buckets,
		uint8_t* types,
		int* summary_suffix_array,
		int* summary_offsets,
		int* final_suffix_array,
		int n
) {
	memset(final_suffix_array, -1, n * sizeof(int));

	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	int  offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	for (int i = n - 1; i >= 0; --i) {
		int str_idx = summary_offsets[summary_suffix_array[i]];

		int char_idx = (int)str[str_idx];
		final_suffix_array[bucket_ends[char_idx]] = str_idx;
		--bucket_ends[char_idx];
	}

	final_suffix_array[0] = n - 1;
}

void construct_sais_suffix_array_int(
	int* str,
	int* suffix_array,
	int n
) {
	memset(suffix_array, -1, n * sizeof(int));

	// STEP 1: Get suffix types
	// Let 0 be S-type and 1 be L-type
	
	// Store n for the empty suffix, which is always S-type.
	uint8_t* types = (uint8_t*)malloc(n);
	types[n - 1] = 0;
	types[n - 2] = 1;
	for (int i = n - 2; i >= 0; --i) {
		if (str[i] == str[i + 1]) {
			types[i] = types[i + 1];
			continue;
		}

		types[i] = (uint8_t)(str[i] > str[i + 1]);
	}

	// STEP 2: Get bucket offsets for chars.
	// We can use counts of each char to find appropriate 
	// starting and ending offsets of all suffixes.
	int* buckets = (int*)malloc(256 * sizeof(int));
	memset(buckets, 0, 256 * sizeof(int));

	for (int i = 0; i < n; ++i) {
		++buckets[str[i]];
	}

	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	memset(bucket_ends, 0, 256 * sizeof(int));

	int* bucket_starts = (int*)malloc(256 * sizeof(int));
	memset(bucket_starts, 0, 256 * sizeof(int));
	// memcpy(bucket_starts, buckets, 256 * sizeof(int));

	// Calculate starting offsets of each char.
	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	// STEP 3: Identify LMS suffixes and use
	// induced sorting to sort them.

	// Iterate through the string to find LMS suffixes.
	// When an LMS suffix is found, place it in the suffix array
	// at the end of its bucket and decrement the bucket end.
	for (int i = 0; i < n - 1; ++i) {
		// If current char is S type and previous is L type, then
		// we have found an LMS char.
		if (types[i] == 0 && types[i - 1] == 1) {
			int bucket_char = str[i];
			int suffix_idx  = bucket_ends[bucket_char];

			suffix_array[suffix_idx] = i;

			// Update the bucket end.
			--bucket_ends[bucket_char];
		}
	}
	suffix_array[0] = n - 1;

	induce_sort_l_chars(str, suffix_array, types, buckets, n);
	induce_sort_s_chars(str, suffix_array, types, buckets, n);

	// STEP 4: Enumerate LMS suffixes.
	int  lms_count    = 0;
	int  lms_offset   = 0;
	int  prev_lms_offset = 0;
	int* lms_suffixes = (int*)malloc(n * sizeof(int));
	memset(lms_suffixes, -1, n * sizeof(int));

	// lms_suffixes[0] = lms_count;
	lms_suffixes[suffix_array[0]] = lms_count;
	prev_lms_offset = suffix_array[0];

	for (int i = 1; i < n; ++i) {
		lms_offset = suffix_array[i];

		// Check if lms suffix
		if (!(types[lms_offset] == 0 && types[lms_offset - 1] == 1)) {
			continue;
		}

		// Compare with previous LMS suffix.
		// Iterate through string at suffix_idx and compare with string 
		// at lms_offset. If they are both equal until the next LMS char
		// is reached, then they are equal and we won't increment lms_count.
		size_t idx = 0;
		bool new_lms = true;
		while (true) {
			bool a_is_lms = (types[lms_offset + idx] == 0 && types[lms_offset + idx - 1] == 1);
			bool b_is_lms = (types[prev_lms_offset + idx] == 0 && types[prev_lms_offset + idx - 1] == 1);

			if (idx > 0 && a_is_lms && b_is_lms) {
				new_lms = false;
				break;
			}

			if (str[lms_offset + idx] != str[prev_lms_offset + idx]) {
				break;
			}

			if (a_is_lms != b_is_lms) {
				break;
			}

			++idx;
		}

		if (new_lms) {
			++lms_count;
		}

		prev_lms_offset 		 = lms_offset;
		lms_suffixes[lms_offset] = lms_count;
	}

	++lms_count;

	// Create summary string.
	int* lms_summary 		 = (int*)malloc(lms_count * sizeof(int));
	int* lms_summary_offsets = (int*)malloc(lms_count * sizeof(int));

	size_t lms_summary_idx = 0;
	for (int i = 0; i < n; ++i) {
		if (lms_suffixes[i] == -1) continue;
		lms_summary[lms_summary_idx] 		 = lms_suffixes[i];
		lms_summary_offsets[lms_summary_idx] = i;
		++lms_summary_idx;
	}

	// Create sorted suffix array for the lms_summary.
	int* lms_summary_suffix_array = (int*)malloc(lms_count * sizeof(int));
	construct_summary_suffix_array(
		lms_summary,
		lms_summary_suffix_array,
		lms_count,
		n
	);

	accurate_lms_sort(
		str,
		buckets,
		types,
		lms_summary_suffix_array,
		lms_summary_offsets,
		suffix_array,
		n
	);

	induce_sort_l_chars(str, suffix_array, types, buckets, n);
	induce_sort_s_chars(str, suffix_array, types, buckets, n);

	print_suffix_array(suffix_array, n);
	return;

	free(types);
	free(buckets);
	free(bucket_ends);
	free(bucket_starts);
	free(lms_suffixes);
}

void construct_summary_suffix_array(
		int* summary_str,
		int* summary_suffix_array,
		int n,
		int summary_alphabet_size
) {
	if (summary_alphabet_size == n - 1) {
		// Bucket sort
		memset(summary_suffix_array, -1, n * sizeof(int));

		summary_suffix_array[0] = n - 1;

		for (int i = 0; i < n - 1; ++i) {
			summary_suffix_array[summary_str[i] + 1] = i;
		}
	}
	else {
		// Recursively construct suffix array for the summary string.
		construct_sais_suffix_array_int(
			summary_str,
			summary_suffix_array,
			n
		);
	}
}



void construct_sais_suffix_array(
	const char* str,
	int* suffix_array,
	int n
) {
	memset(suffix_array, -1, n * sizeof(int));

	// STEP 1: Get suffix types
	// Let 0 be S-type and 1 be L-type
	
	// Store n for the empty suffix, which is always S-type.
	uint8_t* types = (uint8_t*)malloc(n);
	types[n - 1] = 0;
	types[n - 2] = 1;
	for (int i = n - 2; i >= 0; --i) {
		if (str[i] == str[i + 1]) {
			types[i] = types[i + 1];
			continue;
		}

		types[i] = (uint8_t)(str[i] > str[i + 1]);
	}

	// STEP 2: Get bucket offsets for chars.
	// We can use counts of each char to find appropriate 
	// starting and ending offsets of all suffixes.
	int* buckets = (int*)malloc(256 * sizeof(int));
	memset(buckets, 0, 256 * sizeof(int));

	for (int i = 0; i < n; ++i) {
		++buckets[(int)str[i]];
	}

	int* bucket_ends = (int*)malloc(256 * sizeof(int));
	memset(bucket_ends, 0, 256 * sizeof(int));

	int* bucket_starts = (int*)malloc(256 * sizeof(int));
	memset(bucket_starts, 0, 256 * sizeof(int));
	// memcpy(bucket_starts, buckets, 256 * sizeof(int));

	// Calculate starting offsets of each char.
	int offset = 0;
	for (size_t i = 0; i < 256; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
		bucket_ends[i] = offset - 1;
	}

	// STEP 3: Identify LMS suffixes and use
	// induced sorting to sort them.

	// Iterate through the string to find LMS suffixes.
	// When an LMS suffix is found, place it in the suffix array
	// at the end of its bucket and decrement the bucket end.
	for (int i = 0; i < n - 1; ++i) {
		// If current char is S type and previous is L type, then
		// we have found an LMS char.
		if (types[i] == 0 && types[i - 1] == 1) {
			int bucket_char = (int)str[i];
			int suffix_idx  = bucket_ends[bucket_char];

			suffix_array[suffix_idx] = i;

			// Update the bucket end.
			--bucket_ends[bucket_char];
		}
	}
	suffix_array[0] = n - 1;

	induce_sort_l_chars(str, suffix_array, types, buckets, n);
	induce_sort_s_chars(str, suffix_array, types, buckets, n);

	// STEP 4: Enumerate LMS suffixes.
	int  lms_count    = 0;
	int  lms_offset   = 0;
	int  prev_lms_offset = 0;
	int* lms_suffixes = (int*)malloc(n * sizeof(int));
	memset(lms_suffixes, -1, n * sizeof(int));

	// lms_suffixes[0] = lms_count;
	lms_suffixes[suffix_array[0]] = lms_count;
	prev_lms_offset = suffix_array[0];

	for (int i = 1; i < n; ++i) {
		lms_offset = suffix_array[i];

		// Check if lms suffix
		if (!(types[lms_offset] == 0 && types[lms_offset - 1] == 1)) {
			continue;
		}

		// Compare with previous LMS suffix.
		// Iterate through string at suffix_idx and compare with string 
		// at lms_offset. If they are both equal until the next LMS char
		// is reached, then they are equal and we won't increment lms_count.
		size_t idx = 0;
		bool new_lms = true;
		while (true) {
			bool a_is_lms = (types[lms_offset + idx] == 0 && types[lms_offset + idx - 1] == 1);
			bool b_is_lms = (types[prev_lms_offset + idx] == 0 && types[prev_lms_offset + idx - 1] == 1);

			if (idx > 0 && a_is_lms && b_is_lms) {
				new_lms = false;
				break;
			}

			if (str[lms_offset + idx] != str[prev_lms_offset + idx]) {
				break;
			}

			if (a_is_lms != b_is_lms) {
				break;
			}

			++idx;
		}

		if (new_lms) {
			++lms_count;
		}

		prev_lms_offset 		 = lms_offset;
		lms_suffixes[lms_offset] = lms_count;
	}

	++lms_count;

	// Create summary string.
	int* lms_summary 		 = (int*)malloc(lms_count * sizeof(int));
	int* lms_summary_offsets = (int*)malloc(lms_count * sizeof(int));

	size_t lms_summary_idx = 0;
	for (int i = 0; i < n; ++i) {
		if (lms_suffixes[i] == -1) continue;
		lms_summary[lms_summary_idx] 		 = lms_suffixes[i];
		lms_summary_offsets[lms_summary_idx] = i;
		++lms_summary_idx;
	}

	// Create sorted suffix array for the lms_summary.
	int* lms_summary_suffix_array = (int*)malloc(lms_count * sizeof(int));
	construct_summary_suffix_array(
		lms_summary,
		lms_summary_suffix_array,
		lms_count,
		n
	);

	accurate_lms_sort(
		str,
		buckets,
		types,
		lms_summary_suffix_array,
		lms_summary_offsets,
		suffix_array,
		n
	);

	induce_sort_l_chars(str, suffix_array, types, buckets, n);
	induce_sort_s_chars(str, suffix_array, types, buckets, n);

	print_suffix_array(suffix_array, n);
	return;

	free(types);
	free(buckets);
	free(bucket_ends);
	free(bucket_starts);
	free(lms_suffixes);
}
*/


int get_substring_positions(
    char* str,
    uint32_t* suffix_array,
    uint64_t n,
    const char* substring
) {
    int64_t m = strlen(substring);
    int64_t first = 0;
    int64_t last = n - 1;
    int64_t  start = -1, end = -1;

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        int64_t mid = (first + last) / 2;
        if (strncmp(str + suffix_array[mid], substring, m) < 0) {
            first = mid + 1;
        }
        else {
            last = mid - 1;
            start = mid; // Potential start found
        }
    }

    if (start == -1) { // Substring not found
        printf("Substring not found\n");
        return 0;
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
            end = mid; // Potential end found
        }
    }

    // Adjusted to correctly iterate through start to end
	/*
    for (int i = start; i <= end; ++i) {
        printf("Substring found at position %d: %s\n", suffix_array[i], substring);
    }
	*/
	printf("Found %ld matches\n", end - start + 1);

    return end - start + 1;
}


int main() {
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted_1M.csv";
	const char* filename = "/home/jdm365/SearchApp/data/companies_sorted.csv";
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted_100M.csv";

	char* buffer = nullptr;
	uint64_t buffer_size;
	read_text_into_buffer(filename, &buffer, &buffer_size);

	// Just use first 2GB
	// buffer_size = (1 * (uint64_t)1024 * 1024 * 1024) - 1;
	/*
	// char* buffer = (char*)"I work at netflix btw.";
	char* buffer = (char*)"cabbage";
	int buffer_size = strlen(buffer);
	*/

	uint64_t n = buffer_size;
	alignas(64) uint32_t* suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

	auto start = std::chrono::high_resolution_clock::now();
	construct_truncated_suffix_array(buffer, suffix_array, n);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_truncated = end - start;

	/*
	int64_t* suffix_array64 = (int64_t*)malloc(n * sizeof(int64_t));
	int64_t* freq_table = (int64_t*)malloc(256 * sizeof(int64_t));
	start = std::chrono::high_resolution_clock::now();
	libsais64(
		(const uint8_t*)buffer,
		suffix_array64,
		(int64_t)n,
		(int64_t)0,
		freq_table
	);
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	*/

	int* suffix_array2 = (int*)malloc(n * sizeof(int));
	start = std::chrono::high_resolution_clock::now();
	divsufsort(
		(uint8_t*)buffer,
		suffix_array2,
		(int)n
	);
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	printf("Elapsed time construction truncated:  %f seconds\n", elapsed_truncated.count());
	// printf("Elapsed time construction libsais:   %f seconds\n", elapsed.count());
	printf("Elapsed time construction divsufsort: %f seconds\n", elapsed.count());


	uint32_t num_differeces = 0;
	for (uint64_t i = 0; i < n; ++i) {
		if (
				suffix_array[i] == (uint32_t)suffix_array2[i] 
					// ||
				// (suffix_array[i-1] == (uint32_t)suffix_array2[i]) 
					// ||
				// (suffix_array[i-2] == (uint32_t)suffix_array2[i]) 
					// ||
				// (suffix_array[i-3] == (uint32_t)suffix_array2[i]) 
					// ||
				// (suffix_array[i+1] == (uint32_t)suffix_array2[i]) 
			) {
			continue;
		}
		printf("Difference at index %lu: %u, %u\n", i, suffix_array[i], suffix_array2[i]);
		++num_differeces;
	}
	printf("Num differences: %d/%lu\n", num_differeces, n);

	// const char* substr = "jake the god damn snake";
	const char* substr = "ibm";

	start = std::chrono::high_resolution_clock::now();
	int count = get_substring_positions(buffer, suffix_array, n, substr);
	printf("Substring found %d times\n", count);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::micro> elapsed_construction = end - start;

	printf("Elapsed time query: 	   %f microseconds\n", elapsed_construction.count());

	free(buffer);
	free(suffix_array);

	return 0;
}
