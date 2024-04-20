cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native -save-temps -fverbose-asm -g")

project(eisvogel)

option(BUILD_PYTHON "Build with python support." OFF)
option(BUILD_WFC "Also build weighting field creator." OFF)
option(BUILD_TESTS "Build tests." OFF)

add_subdirectory(cmake)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(extern)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)
