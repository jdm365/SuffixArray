# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t, uint64_t
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


cimport numpy as np
import numpy as np
from tqdm import tqdm
import os
from typing import List
np.import_array()

from time import perf_counter


cdef extern from "engine.h":
    void construct_truncated_suffix_array(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
        uint32_t max_suffix_length,
        bool use_index_array
    ) nogil
    void construct_truncated_suffix_array_from_csv(
        const char* csv_file,
        uint32_t column_idx,
        vector[uint32_t]& suffix_array,
        uint32_t* suffix_array_size,
        uint32_t max_suffix_length
    ) nogil
    void _construct_truncated_suffix_array_from_csv(
        const char* csv_file,
        uint32_t column_idx,
        uint32_t* suffix_array,
        uint32_t* suffix_array_size,
        uint32_t max_suffix_length
    ) nogil
    vector[uint32_t] get_matching_indices(
        const char* text,
        uint32_t* suffix_array,
        uint32_t* suffix_array_idxs,
        uint32_t text_length,
        const char* substring,
        int k
    )
    vector[string] get_matching_records(
        const char* filename,
        uint32_t* suffix_array,
        uint32_t text_length,
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
    cdef vector[uint32_t] suffix_array
    cdef vector[uint32_t] suffix_array_idxs
    cdef uint32_t max_suffix_length
    cdef uint32_t text_length
    cdef uint32_t num_rows
    cdef int num_threads
    cdef str csv_file
    cdef str search_column
    cdef uint32_t column_idx 
    cdef str save_dir
    cdef bool use_index_array
    cdef bool from_csv
    cdef list columns

    def __init__(
            self, 
            csv_file: str = '',
            search_column: str = '',
            documents: List[str] = [],
            max_suffix_length: int = 64,
            load_dir: str = '',
            use_index_array: bool = False
            ):
        self.max_suffix_length = <uint32_t>max_suffix_length
        self.use_index_array   = use_index_array
        self.csv_file          = csv_file
        self.search_column     = search_column

        if load_dir != '':
            self.load(load_dir)
            return

        if not isinstance(documents, list):
            try:
                documents = list(documents)
            except:
                raise ValueError("Documents must be a list of strings")

        self.from_csv = False

        if len(documents) > 0:
            self.construct_truncated_suffix_array(documents)
            self.num_rows = len(documents)
        else:
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

            self.from_csv = True
            self.construct_truncated_suffix_array_from_csv()




    def save(self, save_dir: str):
        if save_dir in {'suffix_array', 'tests', 'data'}:
            raise ValueError("Cannot save to reserved directory")

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        else:
            os.system(f"rm -rf {save_dir}/*")

        self.save_dir = save_dir

        cdef FILE* f

        ## Save text
        f = fopen(os.path.join(self.save_dir, 'text.txt').encode('utf-8'), 'wb')
        fwrite(self.text.c_str(), sizeof(char), self.text_length, f)
        fclose(f)

        ## Save suffix array. Just binary dump with c/cpp
        f = fopen(os.path.join(self.save_dir, 'suffix_array.bin').encode('utf-8'), 'wb')
        fwrite(self.suffix_array.data(), sizeof(uint32_t), self.text_length, f)
        fclose(f)

        ## Save suffix array indices
        f = fopen(os.path.join(self.save_dir, 'suffix_array_idxs.bin').encode('utf-8'), 'wb')
        fwrite(self.suffix_array_idxs.data(), sizeof(uint32_t), self.text_length, f)
        fclose(f)

        ## Save metadata
        with open(os.path.join(self.save_dir, 'metadata.txt'), 'w') as file:
            file.write(f"max_suffix_length: {self.max_suffix_length}\n")
            file.write(f"text_length: {self.text_length}\n")
            file.write(f"num_rows: {self.num_rows}\n")
            file.write(f"num_threads: {self.num_threads}\n")


    def load(self, save_dir: str):
        cdef FILE* f

        ## Load text
        f = fopen(os.path.join(save_dir, 'text.txt').encode('utf-8'), 'rb')
        if f is NULL:
            raise ValueError(f"File {save_dir} does not exist")

        fseek(f, 0, SEEK_END)
        cdef uint64_t buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.text.resize(buffer_size)

        fread(self.text.data(), sizeof(char), buffer_size, f)
        fclose(f)

        ## Load suffix array
        f = fopen(os.path.join(save_dir, 'suffix_array.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.suffix_array.resize(buffer_size // sizeof(uint32_t))
        fread(self.suffix_array.data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
        fclose(f)

        ## Load suffix array indices
        f = fopen(os.path.join(save_dir, 'suffix_array_idxs.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)
        self.suffix_array_idxs.resize(buffer_size // sizeof(uint32_t))
        fread(self.suffix_array_idxs.data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
        fclose(f)

        ## Load metadata
        with open(os.path.join(save_dir, 'metadata.txt'), 'r') as file:
            metadata = file.read().split('\n')

            self.max_suffix_length = int(metadata[0].split(': ')[1])
            self.text_length       = int(metadata[1].split(': ')[1])
            self.num_rows          = int(metadata[2].split(': ')[1])
            self.num_threads       = int(metadata[3].split(': ')[1])



    cpdef void construct_truncated_suffix_array(self, list documents):
        ## Convert text to a single string with 
        ## newline separators and get row offsets
        init = perf_counter()
        cdef vector[uint32_t] row_offsets
        self.num_rows = len(documents)

        cdef uint32_t global_offset = 0
        cdef int _id
        row_offsets.push_back(0)
        for _id in range(self.num_rows - 1):
            global_offset += len(documents[_id].encode('utf-8')) + 1
            row_offsets.push_back(global_offset)

        self.text = '\n'.join(documents).encode('utf-8')
        lowercase_string(self.text)
        self.text_length = <uint64_t>len(self.text)

        ## Get the suffix array indices
        self.suffix_array_idxs.resize(self.text_length)
        cdef int i, j
        cdef int n = self.num_rows
        with nogil:
            for i in prange(n-1):
                for j in range(row_offsets[i], row_offsets[i+1]):
                    self.suffix_array_idxs[j] = i

        self.suffix_array_idxs[row_offsets[n-1]] = self.num_rows - 1

        self.suffix_array.resize(self.text_length)
        print(f"Text preprocessed in {perf_counter() - init:.2f} seconds")

        init = perf_counter()
        print("...Constructing suffix array...")
        with nogil:
            construct_truncated_suffix_array(
                self.text.c_str(),
                self.suffix_array.data(),
                self.text_length,
                self.max_suffix_length,
                True
            )
        print(f"Suffix array constructed in {perf_counter() - init:.2f} seconds")


    cpdef void construct_truncated_suffix_array_from_csv(self):

        init = perf_counter()
        print("...Constructing suffix array from csv...")
        cdef string filename = self.csv_file.encode('utf-8')
        with nogil:
            construct_truncated_suffix_array_from_csv(
                filename.c_str(),
                self.column_idx,
                self.suffix_array,
                &self.text_length,
                self.max_suffix_length
            )
        print(f"Suffix array constructed in {perf_counter() - init:.2f} seconds")
        print(f"Text length: {self.text_length}")


    def query(self, substring: str, k: int = 1000):
        cdef np.ndarray[uint32_t, ndim=1] positions = np.array(
            get_matching_indices(
                self.text.c_str(),
                self.suffix_array.data(),
                self.suffix_array_idxs.data(),
                self.text_length,
                substring.lower().encode('utf-8'),
                k
                ),
            dtype=np.uint32
        )
        return positions

    def query_records(self, substring: str, k: int = 1000):
        cdef vector[string] _records = get_matching_records(
            self.csv_file.encode('utf-8'),
            self.suffix_array.data(),
            self.text_length,
            substring.lower().encode('utf-8'),
            k
        )

        records = [x.decode('utf-8').split(',') for x in _records]
        records = [dict(zip(self.columns, x)) for x in records]
        return records
