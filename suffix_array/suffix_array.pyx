# cython: language_level=3

cimport cython

from libc.stdint cimport uint16_t, uint32_t, uint64_t
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
        SEEK_END,

        printf,
        fflush,
        stdout
    )


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
    ## define opaque type
    ctypedef struct buffer_bit
    ctypedef struct SuffixArray:
        uint32_t*   suffix_array
        buffer_bit* is_quoted_bitflag
        uint64_t    global_byte_start_idx
        uint64_t    global_byte_end_idx
        uint32_t    max_suffix_length
        uint32_t    n


    void init_suffix_array(SuffixArray* suffix_array, uint32_t max_suffix_length) nogil
    void init_suffix_array_byte_idxs(
        SuffixArray* suffix_array,
        uint32_t max_suffix_length,
        uint64_t byte_start_idx,
        uint64_t byte_end_idx,
        uint32_t num_bytes
    ) nogil
    void construct_truncated_suffix_array(const char* text, SuffixArray* suffix_array) nogil
    void construct_truncated_suffix_array_from_csv_partitioned(
        const char* csv_file,
        uint32_t column_idx,
        SuffixArray* suffix_array
    ) nogil
    void construct_truncated_suffix_array_from_csv_partitioned_mmap(
        const char* csv_file,
        uint32_t column_idx,
        SuffixArray* suffix_array,
        uint16_t num_columns
    ) nogil
    uint32_t get_matching_records(
        const char* text,
        const SuffixArray* suffix_array,
        const char* substring,
        uint32_t k,
        char** matching_records
    )
    uint32_t get_matching_records_file(
        const char* text,
        const SuffixArray* suffix_array,
        const char* substring,
        uint32_t k,
        char** matching_records
    )

cdef void lowercase_string(char* s, uint64_t length):
    cdef int i
    for i in prange(length, nogil=True):
        if s[i] >= 65 and s[i] <= 90:
            s[i] += 32


cdef class SuffixArrayEngine:
    cdef char*          text
    cdef SuffixArray**  suffix_arrays
    cdef uint32_t       num_rows
    cdef uint32_t       max_suffix_length
    cdef uint32_t       search_col_idx 
    cdef int            num_threads
    cdef str            csv_filename
    cdef str            search_col
    cdef str            save_dir
    cdef list           columns
    cdef uint32_t       num_partitions


    def __init__(self, max_suffix_length: int = 64):
        self.max_suffix_length = max_suffix_length


    cpdef void construct_truncated_suffix_array_documents(self, list documents):
        if not isinstance(documents, list):
            try:
                documents = list(documents)
            except:
                raise ValueError("Documents must be a list of strings")

        self.num_rows = len(documents)
        ## self.partition_byte_boundaries.push_back(0)

        ## Convert text to a single string with 
        ## newline separators and get row offsets
        self.num_rows = len(documents)

        text = '\n'.join(documents).encode('utf-8')
        self.text = text
        lowercase_string(self.text, len(self.text))
        self.text_length = <uint64_t>len(self.text)

        cdef uint32_t TWO_GB    = 2 * 1024 * 1024 * 1024
        cdef uint32_t rem_bytes = self.text_length % TWO_GB

        self.num_partitions = max(1, (self.text_length // TWO_GB) + (rem_bytes > 0))
        self.suffix_arrays  = <SuffixArray**>malloc(self.num_partitions * sizeof(SuffixArray*))

        cdef uint32_t num_bytes
        cdef uint64_t end_byte
        cdef uint64_t byte_idx = 0
        cdef uint64_t idx

        with nogil:
            for idx in range(self.num_partitions):
                num_bytes = TWO_GB if idx < self.num_partitions - 1 else rem_bytes
                end_byte = byte_idx + num_bytes

                init_suffix_array_byte_idxs(
                        self.suffix_arrays[idx], 
                        self.max_suffix_length,
                        byte_idx,
                        end_byte,
                        num_bytes
                        )
                byte_idx += num_bytes

        idx = 0
        while idx < self.num_partitions - 1:
            with nogil:
                construct_truncated_suffix_array(self.text + idx * TWO_GB, self.suffix_arrays[idx])
                idx += 1

        with nogil:
            construct_truncated_suffix_array(self.text + idx * TWO_GB, self.suffix_arrays[idx])


    cpdef void construct_truncated_suffix_array_from_csv(self, filename: str, search_column: str):
        cdef uint64_t idx

        self.csv_filename  = filename
        self.search_col    = search_column

        if not self.csv_filename.endswith('.csv'):
            raise ValueError("Invalid file format. Must be a CSV file")

        ## Check for column name in header.
        ## Get column index.
        with open(self.csv_filename, 'r') as file:
            header = file.readline().strip().split(',')
            if self.search_col not in header:
                raise ValueError(f"Column {self.search_column} not found in CSV file")

            self.columns        = header
            self.search_col_idx = header.index(self.search_col)

        filesize = os.path.getsize(self.csv_filename)
        two_gb   = 2 * 1024 * 1024 * 1024

        self.num_partitions = max(1, (filesize // two_gb) + (filesize % two_gb > 0))
        self.suffix_arrays  = <SuffixArray**>malloc(self.num_partitions * sizeof(SuffixArray*))

        for idx in range(self.num_partitions):
            self.suffix_arrays[idx] = <SuffixArray*>malloc(sizeof(SuffixArray))
            init_suffix_array(self.suffix_arrays[idx], self.max_suffix_length)

        f_name = self.csv_filename.encode('utf-8')
        cdef char* c_filename = f_name

        cdef uint16_t num_columns = len(self.columns)
        with nogil:
            self.suffix_arrays[0].global_byte_start_idx = 0
            for idx in range(self.num_partitions - 1):
                construct_truncated_suffix_array_from_csv_partitioned_mmap(
                    c_filename,
                    self.search_col_idx,
                    self.suffix_arrays[idx],
                    num_columns
                )
                '''
                construct_truncated_suffix_array_from_csv_partitioned(
                    c_filename,
                    self.search_col_idx,
                    self.suffix_arrays[idx]
                )
                '''
                self.suffix_arrays[idx + 1].global_byte_start_idx = self.suffix_arrays[idx].global_byte_end_idx

            '''
            construct_truncated_suffix_array_from_csv_partitioned(
                c_filename,
                self.search_col_idx,
                self.suffix_arrays[idx]
            )

            '''
            construct_truncated_suffix_array_from_csv_partitioned_mmap(
                c_filename,
                self.search_col_idx,
                self.suffix_arrays[idx],
                num_columns
            )



    cpdef query_records(self, substring: str, k: int = 1000):
        if substring == '':
            return []

        cdef uint32_t TWO_GB = 2 * 1024 * 1024 * 1024
        cdef list all_records = []
        cdef char** records = <char**>malloc(k * sizeof(char*))
        cdef int i = 0
        fname = self.csv_filename.encode('utf-8')
        cdef char* c_filename = fname
        cdef uint32_t num_matches = 0

        for i in range(self.num_partitions - 1):

            num_matches = get_matching_records_file(
                    c_filename,
                    self.suffix_arrays[i],
                    substring.lower().encode('utf-8'),
                    k,
                    records
                    )
            for j in range(num_matches):
                all_records.append(records[j].decode('utf-8'))

        num_matches = get_matching_records_file(
                c_filename,
                self.suffix_arrays[self.num_partitions - 1],
                substring.lower().encode('utf-8'),
                k,
                records
                )

        if num_matches == 0:
            free(records)
            return []

        for j in range(num_matches):
            all_records.append(records[j].decode('utf-8'))

        reader      = csv.reader(all_records, delimiter=',')
        all_records = [x for x in reader]
        all_records = [dict(zip(self.columns, x)) for x in all_records]

        for idx in range(num_matches):
            free(records[idx])

        free(records)

        return all_records


    """
    cpdef list query_records(self, substring: str, k: int = 1000):
        cdef list all_records = []
        cdef list records
        cdef vector[string] _records
        cdef int i
        cdef string c_filename = self.csv_file.encode('utf-8')

        print(f"CSV file: {self.csv_file.encode('utf-8')}")
        print(f"Num partitions: {self.num_partitions}")
        print(f"Text lengths: {self.text_lengths}")
        print(f"Partition byte boundaries: {self.partition_byte_boundaries}")
        print(f"Suffix array size: {self.suffix_arrays[0].size()}")
        print(f"Suffix array size: {self.suffix_arrays[1].size()}")
        for i in range(self.num_partitions):
            _records = get_matching_records_file(
                c_filename.c_str(),
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
    """


    """
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
    """
