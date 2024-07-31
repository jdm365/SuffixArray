#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>

// #include "divsufsort.h"
#include "engine.h"
#include "libsais.h"

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


void bench_optimized_suffix_array_construction(
		const char* filename,
		uint32_t max_suffix_length
) {
	SuffixArray suffix_array_struct;
	init_suffix_array(&suffix_array_struct, max_suffix_length);

	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	construct_truncated_suffix_array_from_csv_partitioned(
			filename,
			0,
			&suffix_array_struct
			);
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	double time_taken = (end_time.tv_sec - start_time.tv_sec) + 
						(end_time.tv_nsec - start_time.tv_nsec) / 1e9;

	printf("Elapsed time construction optimized:  %f seconds\n\n\n\n", time_taken);

	free_suffix_array(&suffix_array_struct);
}

void bench_libsais(const char* filename) {
	clock_t start_clock = clock();

	char* buffer = NULL;
	uint64_t buffer_size;
	read_text_into_buffer(filename, &buffer, &buffer_size);

	clock_t end_clock = clock();
	double time_taken_read = (end_clock - start_clock) / (double)CLOCKS_PER_SEC;

	int32_t* suffix_array = (int32_t*)malloc(buffer_size * sizeof(int32_t));
	int32_t* freq_table   = (int32_t*)malloc(256 * sizeof(int32_t));
	int32_t status = libsais(
		(const uint8_t*)buffer,
		suffix_array,
		(int32_t)buffer_size,
		(int32_t)0,
		freq_table
	);
	end_clock = clock();
	double time_taken = (end_clock - start_clock) / (double)CLOCKS_PER_SEC;

	printf("Elapsed time construction libsais:  %f seconds\n\n\n\n", time_taken - time_taken_read);

	free(suffix_array);
	free(freq_table);
	free(buffer);
}



int main() {
	const char* filename = "/home/jdm365/SuffixArray/data/names.csv";
	// const char* filename = "data/english.1024MB";

	const uint32_t max_suffix_length = 32;
	bench_optimized_suffix_array_construction(filename, max_suffix_length);
	bench_libsais(filename);

	return 0;
}
