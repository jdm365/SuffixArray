# cython: language_level=3

cimport cython

from libc.stdint cimport uint32_t 

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp cimport bool

cimport numpy as np
import numpy as np
np.import_array()

from time import perf_counter


cdef extern from "engine.h":



cdef class SuffixArray:
    cdef _SuffixArray* suffix_array 
    cdef list documents


    def __init__(
            self, 
            list documents = []
            ):
        self.documents = documents

    def __cinit__(
            self, 
            *args,
            **kwargs
            ):
        pass
