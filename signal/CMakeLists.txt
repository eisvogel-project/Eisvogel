cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(weighting_field)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -msse2 -ftree-vectorizer-verbose=2 -funroll-loops -march=native -g -fno-omit-frame-pointer -fno-inline -fno-default-inline -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")

add_subdirectory(src)
add_subdirectory(util)
