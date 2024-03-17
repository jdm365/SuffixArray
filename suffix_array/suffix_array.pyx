# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t, uint64_t
from cython.parallel cimport parallel
cimport openmp

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool

cimport numpy as np
import numpy as np
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


cdef class SuffixArray:
    cdef string text
    cdef vector[uint32_t] suffix_array
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

    def __cinit__(
            self, 
            *args,
            **kwargs
            ):
        pass

    def construct_truncated_suffix_array(self, documents):
        ## Convert text to a single string with newline separators
        ## and get row offsets
        cdef vector[uint32_t] row_offsets = [0]
        cdef np.ndarray[uint32_t, ndim=1] str_lengths = np.array(
                [0] + [len(doc.lower().encode('utf-8')) + 1 for doc in documents[:-1]],
                ## [len(doc) + 1 for doc in documents],
                dtype=np.uint32
            )

        self.text = '\n'.join(documents).lower().encode('utf-8')
        self.text_length = <uint64_t>len(self.text)
        self.row_offsets = np.cumsum(str_lengths, dtype=np.uint32).data

        self.suffix_array.resize(self.text_length)

        print("...Constructing suffix array...")
        with nogil:
            construct_truncated_suffix_array(
                self.text.c_str(),
                self.suffix_array.data(),
                self.text_length,
                self.max_suffix_length
            )

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
        '''
        cdef vector[uint32_t] positions = get_matching_indices(
            self.text.c_str(),
            self.suffix_array.data(),
            self.text_length,
            substring.lower().encode('utf-8'),
            self.row_offsets.data(),
            self.num_rows
        )
        '''
        return positions
