set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

if(BUILD_WFC)
    find_package(HDF5 REQUIRED)
    add_executable(testHDF5 testHDF5.cxx)
    target_link_libraries(testHDF5 eisvogel_wfc)
endif(BUILD_WFC)
