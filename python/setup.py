import os, glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("signal",
              sources = glob.glob(os.path.join("${CMAKE_CURRENT_SOURCE_DIR}", "*.pyx")),
              libraries = ["signal"],
              library_dirs = ["${PROJECT_BINARY_DIR}/signal/"],
              language = "c++",
              extra_compile_args = ["-I${CMAKE_CURRENT_SOURCE_DIR}/../src/", "-O3"],
          )
]

for ext_module in ext_modules:
    ext_module.cython_directives = {'language_level': "3"}

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
