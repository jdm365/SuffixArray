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
np.import_array()

import lz4.frame

from time import perf_counter


cdef extern from "engine.h":
    void construct_truncated_suffix_array(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
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
    cdef uint64_t text_length
    cdef uint32_t num_rows
    cdef int num_threads
    cdef str save_dir


    def save(self):
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
        ## f = fopen(os.path.join(self.save_dir, 'suffix_array_idxs.bin').encode('utf-8'), 'wb')
        ## fwrite(self.suffix_array_idxs.data(), sizeof(uint32_t), self.text_length, f)
        ## fclose(f)

        ## Save metadata
        with open(os.path.join(self.save_dir, 'metadata.txt'), 'w') as file:
            file.write(f"max_suffix_length: {self.max_suffix_length}\n")
            file.write(f"text_length: {self.text_length}\n")
            file.write(f"num_rows: {self.num_rows}\n")
            file.write(f"num_threads: {self.num_threads}\n")


    def load(self):
        cdef FILE* f

        ## Load text
        f = fopen(os.path.join(self.save_dir, 'text.txt').encode('utf-8'), 'rb')
        if f is NULL:
            raise ValueError(f"File {self.save_dir} does not exist")

        fseek(f, 0, SEEK_END)
        cdef uint64_t buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.text.resize(buffer_size)

        fread(self.text.data(), sizeof(char), buffer_size, f)
        fclose(f)

        ## Load suffix array
        f = fopen(os.path.join(self.save_dir, 'suffix_array.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)

        self.suffix_array.resize(buffer_size // sizeof(uint32_t))
        fread(self.suffix_array.data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
        fclose(f)

        ## Load suffix array indices
        f = fopen(os.path.join(self.save_dir, 'suffix_array_idxs.bin').encode('utf-8'), 'rb')
        fseek(f, 0, SEEK_END)
        buffer_size = ftell(f)
        fseek(f, 0, SEEK_SET)
        self.suffix_array_idxs.resize(buffer_size // sizeof(uint32_t))
        fread(self.suffix_array_idxs.data(), sizeof(uint32_t), buffer_size // sizeof(uint32_t), f)
        fclose(f)

        ## Load metadata
        with open(os.path.join(self.save_dir, 'metadata.txt'), 'r') as file:
            metadata = file.read().split('\n')

            self.max_suffix_length = int(metadata[0].split(': ')[1])
            self.text_length       = int(metadata[1].split(': ')[1])
            self.num_rows          = int(metadata[2].split(': ')[1])
            self.num_threads       = int(metadata[3].split(': ')[1])


    def __init__(
            self, 
            documents,
            max_suffix_length = 64,
            save_dir: str = 'suffix_array_data'
            ):
        self.max_suffix_length = <uint32_t>max_suffix_length

        ## Create dir and dump binary data into individual files
        if os.path.exists(save_dir):
            raise ValueError(f"File {save_dir} already exists")

        self.save_dir = save_dir
        os.makedirs(self.save_dir)

        if not isinstance(documents, list):
            try:
                documents = list(documents)
            except:
                raise ValueError("Documents must be a list of strings")

        self.construct_truncated_suffix_array(documents)
        self.num_rows = len(documents)


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

        ## Write to disk to avoid holding too much memory.
        cdef FILE* f
        f = fopen(os.path.join(self.save_dir, 'suffix_array_idxs.bin').encode('utf-8'), 'wb')
        fwrite(self.suffix_array_idxs.data(), sizeof(uint32_t), self.text_length, f)
        fclose(f)

        ## Free suffix array indices
        ## self.suffix_array_idxs.clear()
        ## self.suffix_array_idxs.shrink_to_fit()

        self.suffix_array.resize(self.text_length)
        print(f"Text preprocessed in {perf_counter() - init:.2f} seconds")

        init = perf_counter()
        print("...Constructing suffix array...")
        with nogil:
            construct_truncated_suffix_array(
                self.text.c_str(),
                self.suffix_array.data(),
                self.text_length,
                self.max_suffix_length
            )
        print(f"Suffix array constructed in {perf_counter() - init:.2f} seconds")



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
