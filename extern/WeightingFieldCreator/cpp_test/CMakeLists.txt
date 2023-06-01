cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(weighting_field_generator)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(PkgConfig REQUIRED)

pkg_check_modules(MEEP REQUIRED meep)
pkg_check_modules(HDF REQUIRED hdf5)

link_directories(${MEEP_LIBRARY_DIRS} ${HDF_LIBRARY_DIRS})
add_executable(meep_test test.cpp)

target_link_libraries(meep_test ${MEEP_LIBRARIES} ${HDF_LIBRARIES} hdf5_cpp)
target_include_directories(meep_test PUBLIC ${MEEP_INCLUDE_DIRS} ${HDF_INCLUDE_DIRS})
target_compile_options(meep_test PUBLIC ${MEEP_CFLAGS_OTHER} ${HDF_CFLAGS_OTHER})
