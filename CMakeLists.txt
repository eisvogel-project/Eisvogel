cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(eisvogel)

option(BUILD_PYTHON "Build with python support." OFF)
option(BUILD_WFC "Also build weighting field creator." OFF)
option(BUILD_TESTS "Build tests." OFF)

add_subdirectory(cmake)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(extern)

if(BUILD_PYTHON)
    add_subdirectory(cpython)
endif(BUILD_PYTHON)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)
