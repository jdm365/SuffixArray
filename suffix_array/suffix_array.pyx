# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t, uint64_t
from cython.parallel cimport parallel
from cython.parallel cimport prange
cimport openmp

from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool

cimport numpy as np
import numpy as np
from tqdm import tqdm
np.import_array()

from time import perf_counter


cdef extern from "engine.h":
    void construct_truncated_suffix_array(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
        uint32_t max_suffix_length
    ) nogil
    vector[uint32_t] get_substring_positions(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
        const char* substring
    )
    vector[uint32_t] get_matching_indices(
        const char* text,
        uint32_t* suffix_array,
        uint32_t text_length,
        const char* substring,
        const uint32_t* row_offsets,
        uint32_t num_rows,
        int k
    )
    vector[uint32_t] get_matching_indices_fast(
        const char* text,
        uint32_t* suffix_array,
        uint32_t* suffix_array_idxs,
        uint32_t text_length,
        const char* substring,
        int k
    )


cdef class SuffixArray:
    cdef string text
    cdef vector[uint32_t] suffix_array
    cdef vector[uint32_t] suffix_array_idxs
    cdef vector[uint32_t] row_offsets
    cdef uint32_t max_suffix_length
    cdef uint64_t text_length
    cdef uint32_t num_rows
    cdef int num_threads


    def __init__(
            self, 
            documents,
            max_suffix_length = 64
            ):
        self.max_suffix_length = <uint32_t>max_suffix_length
        self.construct_truncated_suffix_array(documents)
        self.num_rows = len(documents)


    cpdef void construct_truncated_suffix_array(self, documents):
        ## Convert text to a single string with newline separators
        ## and get row offsets
        init = perf_counter()
        cdef vector[uint32_t] row_offsets = [0]
        cdef np.ndarray[uint32_t, ndim=1] str_lengths = np.array(
                [0] + [len(doc.lower().encode('utf-8')) + 1 for doc in documents[:-1]],
                dtype=np.uint32
            )
        print(f"Text converted to single string in {perf_counter() - init:.2f} seconds")

        self.text = '\n'.join(documents).lower().encode('utf-8')
        print(f"Text converted to single string in {perf_counter() - init:.2f} seconds")

        self.text_length = <uint64_t>len(self.text)
        self.row_offsets = np.cumsum(str_lengths, dtype=np.uint32).data
        self.num_rows = len(documents)
        print(f"Text converted to single string in {perf_counter() - init:.2f} seconds")

        self.suffix_array.resize(self.text_length)

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

        init = perf_counter()
        ## Get the suffix array indices
        self.suffix_array_idxs.resize(self.text_length)
        cdef int i, j
        cdef int n = self.num_rows
        with nogil:
            for i in prange(n-1):
                for j in range(self.row_offsets[i], self.row_offsets[i+1]):
                    self.suffix_array_idxs[j] = i

        self.suffix_array_idxs[self.row_offsets[n-1]] = self.num_rows - 1
        print(f"Suffix array indices constructed in {perf_counter() - init:.2f} seconds")



    def query(self, substring: str, k: int = 1000):
        cdef np.ndarray[uint32_t, ndim=1] positions = np.array(
            get_matching_indices(
                self.text.c_str(),
                self.suffix_array.data(),
                self.text_length,
                substring.lower().encode('utf-8'),
                self.row_offsets.data(),
                self.num_rows,
                k
                ),
            dtype=np.uint32
        )
        return positions


    def query_idxs(self, substring: str, k: int = 1000):
        cdef np.ndarray[uint32_t, ndim=1] positions = np.array(
            get_matching_indices_fast(
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
