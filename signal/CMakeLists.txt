cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(weighting_field)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -funroll-loops -march=native")

add_subdirectory(src)
add_subdirectory(util)
