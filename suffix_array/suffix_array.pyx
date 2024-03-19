# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t, uint64_t
from cython.parallel cimport prange

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
    ## cdef vector[uint32_t] row_offsets
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

        ## Get the suffix array indices
        self.suffix_array_idxs.resize(self.text_length)
        cdef int i, j
        cdef int n = self.num_rows
        with nogil:
            for i in prange(n-1):
                for j in range(row_offsets[i], row_offsets[i+1]):
                    self.suffix_array_idxs[j] = i

        self.suffix_array_idxs[row_offsets[n-1]] = self.num_rows - 1


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
