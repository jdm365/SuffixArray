#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>


static void print_suffix_array(int* suffix_array, int n) {
	for (int i = 0; i < n; ++i) {
		printf("%d ", suffix_array[i]);
	}
	printf("\n");
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
	int* suffix_array,
	int n
) {
	// Construct simple suffix array from str
	for (int i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

	// Sort suffix array
	std::sort(
		suffix_array,
		suffix_array + n,
		[&](int a, int b) {
			return strcmp(str + a, str + b) < 0;
		}
	);
}

void construct_truncated_suffix_array(
	const char* str,
	int* suffix_array,
	int n,
	int max_suffix_length = 100
) {
	max_suffix_length = std::min(max_suffix_length, n);

	// Construct simple suffix array from str
	for (int i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}

	// Sort suffix array only up to max_suffix_length
	std::sort(
		suffix_array,
		suffix_array + n,
		[&](int a, int b) {
			return strncmp(str + a, str + b, max_suffix_length) < 0;
		}
	);

	for (int i = 0; i < n; ++i) {
		printf("%d ", suffix_array[i]);
	}
}

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

	/*
	// Print LMS suffixes
	for (int i = 0; i < lms_count; ++i) {
		printf("%d ", (int)lms_summary[i]);
	}
	printf("\n");
	*/


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

	/*
	// Print LMS suffixes
	for (int i = 0; i < lms_count; ++i) {
		printf("%d ", (int)lms_summary[i]);
	}
	printf("\n");
	*/


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


int get_substring_positions(
    char* str,
    int* suffix_array,
    int n,
    const char* substring
) {
    int m = strlen(substring);
    int first = 0;
    int last = n - 1;
    int start = -1, end = -1;

    // Binary search for the first occurrence of the substring
    while (first <= last) {
        int mid = (first + last) / 2;
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
        int mid = (first + last) / 2;
        if (strncmp(str + suffix_array[mid], substring, m) > 0) {
            last = mid - 1;
        }
        else {
            first = mid + 1;
            end = mid; // Potential end found
        }
    }

    // Adjusted to correctly iterate through start to end
    for (int i = start; i <= end; ++i) {
        // printf("Substring found at position %d: %s\n", suffix_array[i], str + suffix_array[i]);
        printf("Substring found at position %d: %s\n", suffix_array[i], substring);
    }

    return end - start + 1;
}


int main() {
	const char* filename = "/home/jdm365/SearchApp/basic_search/data/companies_sorted_100k.csv";
	// const char* filename = "/home/jdm365/SearchApp/basic_search/data/companies_sorted.csv";

	/*
	char* buffer = nullptr;
	uint64_t buffer_size;
	read_text_into_buffer(filename, &buffer, &buffer_size);
	*/
	// char* buffer = (char*)"I work at netflix btw.";
	char* buffer = (char*)"cabbage";
	int buffer_size = strlen(buffer);

	int n = buffer_size + 1;
	int* suffix_array = (int*)malloc(n * sizeof(int));

	auto start = std::chrono::high_resolution_clock::now();
	// construct_suffix_array(buffer, suffix_array, n);
	construct_truncated_suffix_array(buffer, suffix_array, n);
	printf("\n\n");
	construct_sais_suffix_array(buffer, suffix_array, n);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	const char* substr = "age";
	// const char* substr = "netflix";

	start = std::chrono::high_resolution_clock::now();
	int count = get_substring_positions(buffer, suffix_array, n, substr);
	printf("Substring found %d times\n", count);

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::micro> elapsed_construction = end - start;

	printf("Elapsed time construction: %f seconds\n", elapsed.count());
	printf("Elapsed time query: 	   %f microseconds\n", elapsed_construction.count());

	// free(buffer);
	free(suffix_array);

	return 0;
}
