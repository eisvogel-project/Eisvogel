import os, glob
from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("eisvogel",
              sources = [os.path.join("${CMAKE_CURRENT_SOURCE_DIR}", "eisvogel.pyx")],
              libraries = ["eisvogel"],
              library_dirs = ["${PROJECT_BINARY_DIR}/eisvogel/"],
              language = "c++",
              extra_compile_args = ["-I${PROJECT_SOURCE_DIR}/include/", "-I${PROJECT_SOURCE_DIR}/src/", "-I${MPI_INCLUDE_PATH}",
                                    "-O3", "-std=gnu++2a",
                                    "-ftree-vectorize", "-ffast-math",
                                    "-ftree-vectorizer-verbose=2", "-funroll-loops", "-march=native"],
          )
]

extensions = cythonize(
    extensions,
    compiler_directives = {
        "language_level": "3"
    }
)

setup(
    name = "eisvogel",
    ext_modules = extensions
)
