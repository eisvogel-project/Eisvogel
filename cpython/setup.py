import os, glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("eisvogel",
              sources = [os.path.join("${CMAKE_CURRENT_SOURCE_DIR}", "eisvogel.pyx")],
              libraries = ["eisvogel"],
              library_dirs = ["${PROJECT_BINARY_DIR}/eisvogel/"],
              language = "c++",
              extra_compile_args = ["-I${PROJECT_SOURCE_DIR}/include/", "-O3", "-std=gnu++2a",
                                    "-ftree-vectorize", "-ffast-math", "-msse2", 
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
