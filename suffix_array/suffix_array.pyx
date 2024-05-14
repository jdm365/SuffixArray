# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t, uint64_t
from libc.stdio cimport printf
from cython.parallel cimport prange
from libc.stdlib cimport malloc, free
from libc.stdio cimport (
        fopen, 
        fclose, 
        fwrite, 
        fread, 
        fseek, 
        ftell, 
        FILE, 
        SEEK_SET, 
        SEEK_END
    )

from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool

from cpython.bytes cimport PyBytes_AsString


cimport numpy as np
import numpy as np
from tqdm import tqdm
import os
import sys
import csv
from typing import List
np.import_array()

from time import perf_counter


cdef extern from "engine.h":
    void construct_truncated_suffix_array(
        const char* text,
        vector[uint32_t]& suffix_array,
        uint32_t text_length,
        uint32_t max_suffix_length
    ) nogil
    void construct_truncated_suffix_array_from_csv_partitioned(
        const char* csv_file,
        uint32_t column_idx,
        vector[uint32_t]& suffix_array,
        uint32_t* suffix_array_size,
        uint32_t max_suffix_length,
        uint64_t start_idx,
        uint64_t& end_idx
    ) nogil
    vector[string] get_matching_records(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
        const char* substring,
        int k
    )
    vector[string] get_matching_records_file(
        const char* filename,
        uint32_t* suffix_array,
        uint32_t text_length,
        uint64_t start_idx,
        const char* substring,
        int k
    )

cdef void lowercase_string(string& s) nogil:
    cdef int i
    for i in prange(s.size()):
        if s[i] >= 65 and s[i] <= 90:
            s[i] += 32



cdef class SuffixArray:
    cdef string text
    cdef vector[vector[uint32_t]] suffix_arrays
    cdef uint64_t text_length
    cdef vector[uint32_t] text_lengths
    cdef uint32_t num_rows
    cdef uint32_t max_suffix_length
    cdef int num_threads
    cdef str csv_file
    cdef str search_column
    cdef uint32_t column_idx 
    cdef str save_dir
    cdef bool from_csv
    cdef list columns
    cdef int  num_partitions
    cdef vector[uint64_t] partition_byte_boundaries


    def __init__(self, max_suffix_length: int = 64):
        self.max_suffix_length = max_suffix_length


    cpdef void construct_truncated_suffix_array_documents(self, list documents):
        if not isinstance(documents, list):
            try:
                documents = list(documents)
            except:
                raise ValueError("Documents must be a list of strings")

        self.num_rows = len(documents)
        self.num_partitions = 1
        self.partition_byte_boundaries.push_back(0)

        ## Convert text to a single string with 
        ## newline separators and get row offsets
        self.num_rows = len(documents)

        self.text = '\n'.join(documents).encode('utf-8')
        lowercase_string(self.text)
        self.text_length = <uint64_t>len(self.text)

        cdef uint32_t TWO_GB = 2 * 1024 * 1024 * 1024
        self.num_partitions = max(1, (self.text_length // TWO_GB) + (self.text_length % TWO_GB > 0))
        cdef uint32_t rem_bytes = self.text_length % TWO_GB
        self.suffix_arrays.resize(self.num_partitions)

        cdef uint32_t i = 0
        while i < self.num_partitions - 1:
            with nogil:
                self.suffix_arrays[i].resize(TWO_GB)
                construct_truncated_suffix_array(
                    self.text.data() + i * TWO_GB,
                    self.suffix_arrays[i],
                    TWO_GB,
                    self.max_suffix_length
                )
                i += 1

        with nogil:
            self.suffix_arrays[i].resize(rem_bytes)
            construct_truncated_suffix_array(
                self.text.data() + i * TWO_GB,
                self.suffix_arrays[i],
                rem_bytes,
                self.max_suffix_length
            )


    cpdef void construct_truncated_suffix_array_from_csv(self, filename: str, search_column: str):
        self.csv_file          = filename
        self.search_column     = search_column
        self.from_csv          = True

        if not self.csv_file.endswith('.csv'):
            raise ValueError("Invalid file format. Must be a CSV file")

        ## Check for column name in header.
        ## Get column index.
        with open(self.csv_file, 'r') as file:
            header = file.readline().strip().split(',')
            if self.search_column not in header:
                raise ValueError(f"Column {self.search_column} not found in CSV file")

            self.columns = header
            self.column_idx = header.index(self.search_column)

        filesize = os.path.getsize(self.csv_file)
        two_gb = 2 * 1024 * 1024 * 1024

        self.num_partitions = max(1, (filesize // two_gb) + (filesize % two_gb > 0))

        self.suffix_arrays.resize(self.num_partitions)
        self.text_lengths.resize(self.num_partitions)
        self.partition_byte_boundaries.resize(self.num_partitions + 1)
        self.partition_byte_boundaries[0] = 0


        cdef const char* c_filename = PyBytes_AsString(self.csv_file.encode('utf-8'))
        cdef int i

        with nogil:
            for i in range(self.num_partitions):
                construct_truncated_suffix_array_from_csv_partitioned(
                    c_filename,
                    self.column_idx,
                    self.suffix_arrays[i],
                    &self.text_lengths[i],
                    self.max_suffix_length,
                    self.partition_byte_boundaries[i],
                    self.partition_byte_boundaries[i + 1]
                )



    cpdef query_records_2(self, substring: str, k: int = 1000):
        cdef uint32_t TWO_GB = 2 * 1024 * 1024 * 1024
        cdef list all_records = []
        cdef vector[string] records
        cdef int i = 0

        while i < self.num_partitions - 1:

            records = get_matching_records(
                    self.text.data() + i * TWO_GB,
                    self.suffix_arrays[i].data(),
                    TWO_GB,
                    substring.lower().encode('utf-8'),
                    k
                    )
            all_records.extend(records)
            i += 1

        records = get_matching_records(
                self.text.data() + i * TWO_GB,
                self.suffix_arrays[i].data(),
                self.text_length % TWO_GB,
                substring.lower().encode('utf-8'),
                k
                )
        all_records.extend(records)
        return [x.decode('utf-8') for x in all_records]


    cpdef list query_records(self, substring: str, k: int = 1000):
        cdef list all_records = []
        cdef list records
        cdef vector[string] _records
        cdef int i

        print(f"CSV file: {self.csv_file.encode('utf-8')}")
        print(f"Num partitions: {self.num_partitions}")
        print(f"Text lengths: {self.text_lengths}")
        print(f"Partition byte boundaries: {self.partition_byte_boundaries}")
        print(f"Suffix array size: {self.suffix_arrays[0].size()}")
        print(f"Suffix array size: {self.suffix_arrays[1].size()}")
        for i in range(self.num_partitions):
            _records = get_matching_records_file(
                self.csv_file.encode('utf-8'),
                self.suffix_arrays[i].data(),
                self.text_lengths[i],
                self.partition_byte_boundaries[i],
                substring.lower().encode('utf-8'),
                k
            )

            records = [x.decode('utf-8') for x in _records]

            reader = csv.reader(records, delimiter=',')
            records = [x for x in reader]
            records = [dict(zip(self.columns, x)) for x in records]

            all_records.extend(records)

            if len(all_records) >= k:
                all_records = all_records[:k]
                break

        return all_records


    def save(self, save_dir: str):
        if save_dir in {'suffix_array', 'tests', 'data'}:
            raise ValueError("Cannot save to reserved directory")

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        else:
            os.system(f"rm -rf {save_dir}/*")

        self.save_dir = save_dir

        cdef FILE* f

        ## Save suffix arrays. Just binary dump with c/cpp
        for i in range(self.num_partitions):
            f = fopen(os.path.join(self.save_dir, f'suffix_array_{i}.bin').encode('utf-8'), 'wb')
            fwrite(self.suffix_arrays[i].data(), sizeof(uint32_t), self.text_lengths[i], f)
            fclose(f)

        ## Save text_lengths and partition_byte_boundaries vectors
        f = fopen(os.path.join(self.save_dir, 'text_lengths.bin').encode('utf-8'), 'wb')
        fwrite(self.text_lengths.data(), sizeof(uint32_t), self.num_partitions, f)
        fclose(f)

        f = fopen(os.path.join(self.save_dir, 'partition_byte_boundaries.bin').encode('utf-8'), 'wb')
        fwrite(self.partition_byte_boundaries.data(), sizeof(uint64_t), self.num_partitions + 1, f)
        fclose(f)

        ## Save metadata
        with open(os.path.join(self.save_dir, 'metadata.txt'), 'w') as file:
            file.write(f"max_suffix_length: {self.max_suffix_length}\n")
            file.write(f"text_length: {self.text_length}\n")
            file.write(f"num_rows: {self.num_rows}\n")
            file.write(f"num_threads: {self.num_threads}\n")
            file.write(f"num_partitions: {self.num_partitions}\n")
            if self.from_csv:
                file.write(f"csv_file: {self.csv_file}\n")
            else:
                file.write(f"from_csv: {0}\n")

        ## Save text
        if not self.from_csv:
            f = fopen(os.path.join(self.save_dir, 'text.txt').encode('utf-8'), 'wb')
            fwrite(self.text.c_str(), sizeof(char), self.text_length, f)
            fclose(f)


    def load(self, save_dir: str):
        ## Load metadata
        with open(os.path.join(save_dir, 'metadata.txt'), 'r') as file:
            metadata = file.read().split('\n')

            self.max_suffix_length = int(metadata[0].split(': ')[1])
            self.text_length       = int(metadata[1].split(': ')[1])
            self.num_rows          = int(metadata[2].split(': ')[1])
            self.num_threads       = int(metadata[3].split(': ')[1])
            self.num_partitions    = int(metadata[4].split(': ')[1])
            self.csv_file          = metadata[5].split(': ')[1]
            self.from_csv          = self.csv_file != '0'
        
        print(f"Loading from CSV: {self.from_csv}")
        print(f"CSV file: {self.csv_file}")
        sys.stdout.flush()

        cdef FILE* f
        cdef uint64_t buffer_size

        if not self.from_csv:
            ## Load text
            f = fopen(os.path.join(save_dir, 'text.txt').encode('utf-8'), 'rb')
            if f is NULL:
                raise ValueError(f"File {save_dir} does not exist")

            fseek(f, 0, SEEK_END)
            buffer_size = ftell(f)
            fseek(f, 0, SEEK_SET)

            self.text.resize(buffer_size)

            fread(self.text.data(), sizeof(char), buffer_size, f)
            fclose(f)

        ## Load suffix arrays
        self.suffix_arrays.resize(self.num_partitions)
        for i in range(self.num_partitions):
            f = fopen(os.path.join(save_dir, f'suffix_array_{i}.bin').encode('utf-8'), 'rb')
            fseek(f, 0, SEEK_END)
            buffer_size = ftell(f)
            fseek(f, 0, SEEK_SET)

            self.suffix_arrays[i].resize(buffer_size // sizeof(uint32_t))
            fread(self.suffix_arrays[i].data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
            fclose(f)

        ## Load text_lengths and partition_byte_boundaries vectors
        f = fopen(os.path.join(save_dir, 'text_lengths.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.text_lengths.resize(buffer_size // sizeof(uint32_t))
        fread(self.text_lengths.data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
        fclose(f)

        f = fopen(os.path.join(save_dir, 'partition_byte_boundaries.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.partition_byte_boundaries.resize(buffer_size // sizeof(uint64_t))
        fread(self.partition_byte_boundaries.data(), sizeof(uint64_t), buffer_size // sizeof(uint64_t), f)
        fclose(f)
