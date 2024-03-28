set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native -save-temps -fverbose-asm -g")

add_executable(testSparseNDArray testSparseNDArray.cxx)
target_link_libraries(testSparseNDArray eisvogel)

add_executable(testDenseNDArray testDenseNDArray.cxx)
target_link_libraries(testDenseNDArray eisvogel)

add_executable(testArrayOps testArrayOps.cxx)
target_link_libraries(testArrayOps eisvogel)

add_executable(testVector testVector.cxx)
target_link_libraries(testVector eisvogel)

add_executable(testCompression testCompression.cxx)
target_link_libraries(testCompression eisvogel)

add_executable(testDistributedNDVecArray testDistributedNDVecArray.cxx)
target_link_libraries(testDistributedNDVecArray eisvogel)

add_executable(testCache testCache.cxx)
target_link_libraries(testCache eisvogel)
