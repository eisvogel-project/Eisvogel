cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native -save-temps -fverbose-asm -g")

project(eisvogel)

option(BUILD_TESTS "Build tests." OFF)
option(BUILD_MEEP "Also build weighting field creator." OFF)
option(BUILD_EXTERN "Build external tools." OFF)

add_subdirectory(eisvogel)
add_subdirectory(examples)

if(BUILD_EXTERN)
    add_subdirectory(extern)
endif(BUILD_EXTERN)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)
