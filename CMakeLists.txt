cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(eisvogel)

option(BUILD_PYTHON "Build with python support." OFF)
option(BUILD_WFC "Also build weighting field creator." OFF)

add_subdirectory(cmake)
add_subdirectory(src)
if(BUILD_PYTHON)
    add_subdirectory(cpython)
endif(BUILD_PYTHON)
add_subdirectory(examples)
