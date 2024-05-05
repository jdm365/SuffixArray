#include <cstddef>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <stdbool.h>
#include <omp.h>

#include <iostream>
#include <string>

// #include "divsufsort.h"
#include "engine.h"

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


/*
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

	size_t _ = fread(*buffer, 1, *buffer_size, file);
	fclose(file);

	fflush(stdout);
}
*/


void construct_suffix_array_default(
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


/*
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
*/


void custom_strncmp(const char* str1, const char* str2, int max_len, int& result) {
    for (int i = 0; i < max_len; ++i) {
        if (str1[i] == '\n' || str2[i] == '\n' || str1[i] != str2[i]) {
            result = (unsigned char)str1[i] - (unsigned char)str2[i];
            return;
        }
    }
    result = 0;
}

void recursive_bucket_sort_index(
		const char* str,
		uint32_t* suffix_array,
		uint32_t* temp_suffix_array,
		int string_length,
		uint64_t n,
		int current_depth
) {
	if (n < 32) {
		// Do insertion sort
		for (uint64_t i = 1; i < n; ++i) {
            int j = i;
            while (j > 0) {
                int cmp_result;
                custom_strncmp(
                    str + suffix_array[j] + current_depth,
                    str + suffix_array[j - 1] + current_depth,
                    string_length - (suffix_array[j] + current_depth),
                    cmp_result
                );
                if (cmp_result < 0) {
                    std::swap(suffix_array[j], suffix_array[j - 1]);
                }
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

	for (uint32_t i = 0; i < n; ++i) {
		int char_idx = suffix_array[i] + current_depth;
		uint8_t char_val = char_idx < string_length && str[char_idx] != '\n' ? (uint8_t)str[char_idx] : 0;
		++buckets[char_val];
	}

	uint64_t offset = 0;
	for (size_t i = 0; i < NUM_BUCKETS; ++i) {
		bucket_starts[i] = offset;
		offset += buckets[i];
	}


	memcpy(temp_suffix_array, suffix_array, n * sizeof(uint32_t));
    for (uint64_t i = 0; i < n; ++i) {
        int char_idx = suffix_array[i] + current_depth;
        if (char_idx < string_length && str[char_idx] != '\n') {
            uint8_t char_val = (uint8_t)str[char_idx];
            temp_suffix_array[bucket_starts[char_val]] = suffix_array[i];
            ++bucket_starts[char_val];
        } else {
            // Treat newline or end of string separately if needed
            temp_suffix_array[bucket_starts[0]] = suffix_array[i];
            ++bucket_starts[0];
        }
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

			recursive_bucket_sort_index(
				str,
				suffix_array + bucket_start,
				temp_suffix_array + bucket_start,
				string_length,
				bucket_end - bucket_start,
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

		recursive_bucket_sort_index(
			str,
			suffix_array + bucket_start,
			temp_suffix_array + bucket_start,
			string_length,
			bucket_end - bucket_start,
			current_depth + 1
		);
	}
}

/*
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
}
*/

void construct_suffix_array_index(
	const char* str,
	uint32_t* suffix_array,
	uint64_t n,
	const char delim_char = '\n'
) {
	// Construct simple suffix array from str
	for (uint64_t i = 0; i < n; ++i) {
		suffix_array[i] = i;
	}
	alignas(64) uint32_t* temp_suffix_array = (uint32_t*)malloc(n * sizeof(uint32_t));

	recursive_bucket_sort_index(
		str,
		suffix_array,
		temp_suffix_array,
		n,
		n,
		0
	);
	free(temp_suffix_array);
}


/*
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
	// for (int j = 0; j < 9; ++j) {
		// printf("SUBSTRING: ");
		// for (int i = 0; i < 25; ++i) {
			// if (str[suffix_array[start + j] + i] == '\n') {
				// break;
			// }
			// printf("%c", str[suffix_array[start + j] + i]);
		// }
		// printf("\n\n");
	// }
	// printf("Found %ld matches\n", end - start + 1);

    return end - start + 1;
}
*/


int main() {
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted_100k.csv";
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted_1M.csv";
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted.csv";
	// const char* filename = "/home/jdm365/SearchApp/data/companies_sorted_100M.csv";
	const char* filename = "/home/jdm365/SearchApp/data/names_only.txt";
	// const char* filename = "/home/jdm365/SearchApp/data/names_only_100M.txt";
	// const char* filename = "data/english.1024MB";

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

	uint32_t max_suffix_length = 32;

	/*
	auto start = std::chrono::high_resolution_clock::now();
	construct_truncated_suffix_array(
			buffer, 
			suffix_array, 
			n, 
			max_suffix_length,
			false
			);
	// construct_suffix_array_index(buffer, suffix_array, n);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_truncated = end - start;

	printf("Elapsed time construction truncated:  %f seconds\n\n\n\n", elapsed_truncated.count());
	fflush(stdout);
	*/

	const char* FILENAME = "/home/jdm365/SearchApp/data/companies_sorted.csv";

	auto start = std::chrono::high_resolution_clock::now();
	// start = std::chrono::high_resolution_clock::now();
	uint32_t suffix_array_size;
	_construct_truncated_suffix_array_from_csv(
			FILENAME,
			0, 
			&suffix_array, 
			&suffix_array_size,
			max_suffix_length
			);
	// construct_suffix_array_index(buffer, suffix_array, n);
	auto end = std::chrono::high_resolution_clock::now();
	// end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_truncated = end - start;
	// elapsed_truncated = end - start;

	printf("Elapsed time construction truncated:  %f seconds\n\n\n\n", elapsed_truncated.count());
	fflush(stdout);

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

	/*
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
	*/


	/*
	uint32_t num_differeces = 0;
	for (uint64_t i = 0; i < n; ++i) {
		if (suffix_array[i] == (uint32_t)suffix_array2[i]) {
			continue;
		}
		printf("Difference at index %lu: %u, %u\n", i, suffix_array[i], suffix_array2[i]);
		++num_differeces;
	}
	printf("Num differences: %d/%lu\n", num_differeces, n);
	*/

	// const char* substr = "JAKE THE SNAKE";
	const char* substr = "netflix";
	// const char* substr = "THE";
	/*
	FILE* file = fopen(FILENAME, "r");
	if (file == NULL) {
		printf("Error: File not found\n");
		exit(1);
	}
	char* buffer2 = NULL;
	uint64_t buffer_size2;
	read_text_into_buffer(FILENAME, &buffer2, &buffer_size2);
	n = buffer_size2;

	start = std::chrono::high_resolution_clock::now();
	// std::pair<uint32_t, uint32_t> range = get_substring_positions(buffer, suffix_array, n, substr);
	std::vector<uint32_t> idxs = get_matching_indices_no_idxs(buffer2, suffix_array, n, substr, 10);
	for (uint32_t i = 0; i < idxs.size(); ++i) {
		printf("Substring found at index %u: ", idxs[i]);
		for (uint32_t j = 0; j < 250; ++j) {
			if (buffer2[idxs[i] + j] == '\n') {
				break;
			}
			printf("%c", buffer2[idxs[i] + j]);
		}
		printf("\n");
		fflush(stdout);
	}
	*/
	start = std::chrono::high_resolution_clock::now();

	std::vector<std::string> records = get_matching_records(
			FILENAME,
			suffix_array, 
			(uint64_t)suffix_array_size,
			substr, 
			10
			);
	for (uint32_t i = 0; i < records.size(); ++i) {
		std::cout << records[i] << std::endl;
	}

	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::micro> elapsed_construction = end - start;

	printf("Elapsed time query: 	   %f microseconds\n", elapsed_construction.count());

	free(buffer);
	free(suffix_array);

	return 0;
}
