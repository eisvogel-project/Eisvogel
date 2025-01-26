cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native -save-temps -fverbose-asm -g")

project(eisvogel)

option(BUILD_EXAMPLES "Build examples." OFF)
option(BUILD_TESTS "Build tests." OFF)
option(BUILD_MEEP "Also build weighting field creator." OFF)
option(BUILD_CORSIKA "Build with C8 support." OFF)
option(BUILD_EXTERN "Build external tools." OFF)

add_subdirectory(eisvogel)
add_subdirectory(applications)

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)

if(BUILD_EXTERN OR BUILD_CORSIKA)
  add_subdirectory(extern)
endif(BUILD_EXTERN OR BUILD_CORSIKA)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(BUILD_TESTS)
