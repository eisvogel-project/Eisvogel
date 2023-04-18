import os, glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("pyeisvogel.kernels",
              sources = [os.path.join("${CMAKE_CURRENT_SOURCE_DIR}", "kernels.pyx")],
              libraries = ["eisvogel"],
              library_dirs = ["${PROJECT_BINARY_DIR}/eisvogel/"],
              language = "c++",
              extra_compile_args = ["-I${PROJECT_SOURCE_DIR}/include/", "-O3"],
          ),
    Extension("pyeisvogel.wf",
              sources = [os.path.join("${CMAKE_CURRENT_SOURCE_DIR}", "weightingfield.pyx")],
              libraries = ["eisvogel"],
              library_dirs = ["${PROJECT_BINARY_DIR}/eisvogel/"],
              language = "c++",
              extra_compile_args = ["-I${PROJECT_SOURCE_DIR}/include/", "-O3", "-std=c++20"],
          )
]

extensions = cythonize(
    extensions,
    compiler_directives = {
        "language_level": "3"
    }
)

setup(
    name = "pyeisvogel",
    ext_modules = extensions
)
