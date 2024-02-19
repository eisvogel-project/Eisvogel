set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

add_executable(testSparseNDArray testSparseNDArray.cxx)
target_link_libraries(testSparseNDArray eisvogel)

add_executable(testDenseNDArray testDenseNDArray.cxx)
target_link_libraries(testDenseNDArray eisvogel)

