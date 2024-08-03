from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

MODULE_NAME = "suffix_array"

COMPILER_FLAGS = [
    ## "-std=c++17",
    ## "-stdlib=libc++",
    "-O3",
    "-march=native",
    "-ffast-math",
    "-Wno-unused-variable",
    "-Wno-unknown-pragmas",
    "-Wno-unused-result",
    "-fopenmp",
    ## "-fPIC",
]

LINK_ARGS = [
    ## "-lc++",
    ## "-lc++abi",
    "-L/usr/local/lib",
    "-lgomp",
    ]


extensions = [
    Extension(
        MODULE_NAME,
        sources=["suffix_array/suffix_array.pyx", "suffix_array/engine.c"],
        extra_compile_args=COMPILER_FLAGS,
        language="c",
        include_dirs=[np.get_include(), "suffix_array"],
        extra_link_args=LINK_ARGS
    ),
]

setup(
    name=MODULE_NAME,
    ext_modules=cythonize(extensions),
)
