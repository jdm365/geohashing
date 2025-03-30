from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import os

os.environ["CC"] = "/opt/homebrew/opt/llvm/bin/clang"

MODULE_NAME = "geohashing"

COMPILER_FLAGS = [
        "-O3",
        "-Wall",
        "-march=native",
        "-ffast-math",
        "-isysroot", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
        ## "-I/opt/homebrew/opt/libomp/include",
        ## "-fopenmp",
        "-pthread",
        ]
LINKER_FLAGS = [
    "-L/Library/Developer/CommandLineTools/usr/lib",
    "-isysroot", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
    ## "-L/opt/homebrew/opt/libomp/lib",
    ## "-lomp",
]

extensions = [
        Extension(
            MODULE_NAME,
            sources=[f"{MODULE_NAME}/*.pyx", f"{MODULE_NAME}/engine.c"],
            include_dirs=[
                MODULE_NAME, 
                np.get_include(),
                ],
            extra_compile_args=COMPILER_FLAGS,
            extra_link_args=LINKER_FLAGS,
            language="c",
            ),
        ]

setup(name=MODULE_NAME, ext_modules=cythonize(extensions))
