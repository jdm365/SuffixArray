from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import os

MODULE_NAME = "suffix_array"

COMPILER_FLAGS = [
    "-std=c++17",
    ## "-stdlib=libc++",
    "-O3",
    "-Wall",
    "-Wextra",
    "-march=native",
    "-ffast-math",
    ## "-fopenmp",
    ## "-fPIC",
]

LINK_ARGS = [
    "-lc++",
    "-lc++abi",
    "-L/usr/local/lib"
    ]

OS = os.uname().sysname
CXX = os.environ.get('CXX')

if OS != 'Darwin':
    COMPILER_FLAGS += ['-fopenmp']
    LINK_ARGS      += ['-fopenmp']

if CXX is None:
    raise ValueError("CXX environment variable is not set")
elif CXX == "clang++":
    COMPILER_FLAGS += ["-stdlib=libc++"]
    LINK_ARGS      += ["-stdlib=libc++"]

extensions = [
    Extension(
        MODULE_NAME,
        sources=["suffix_array/suffix_array.pyx", "suffix_array/engine.cpp"],
        extra_compile_args=COMPILER_FLAGS,
        language="c++",
        include_dirs=[np.get_include(), "suffix_array"],
        extra_link_args=LINK_ARGS
    ),
]

setup(
    name=MODULE_NAME,
    ext_modules=cythonize(extensions),
)
