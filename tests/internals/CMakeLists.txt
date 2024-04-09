set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native -save-temps -fverbose-asm -g")
# set(CMAKE_CXX_FLAGS "-O2 -Wall -Wpedantic -Wextra -ftree-vectorize -ffast-math -fno-trapping-math -ftree-vectorizer-verbose=2 -march=native -save-temps -fverbose-asm -g -DNDEBUG")

# add_executable(testSparseNDArray testSparseNDArray.cxx)
# target_link_libraries(testSparseNDArray eisvogel)

# add_executable(testDenseNDArray testDenseNDArray.cxx)
# target_link_libraries(testDenseNDArray eisvogel)

# add_executable(testArrayOps testArrayOps.cxx)
# target_link_libraries(testArrayOps eisvogel)

# add_executable(testVector testVector.cxx)
# target_link_libraries(testVector eisvogel)

# add_executable(testNDVecArray testNDVecArray.cxx)
# target_link_libraries(testNDVecArray eisvogel)

# add_executable(testIndexing testIndexing.cxx)
# target_link_libraries(testIndexing eisvogel)

# add_executable(testCache testCache.cxx)
# target_link_libraries(testCache eisvogel)

# add_executable(testInterpolation testInterpolation.cxx)

add_executable(testCompression testCompression.cxx)
target_link_libraries(testCompression eisvogel)

add_executable(testDistributedNDVecArray testDistributedNDVecArray.cxx)
target_link_libraries(testDistributedNDVecArray eisvogel)

# add_executable(testNDVecArrayInterpolation testNDVecArrayInterpolation.cxx)
# target_link_libraries(testNDVecArrayInterpolation eisvogel)

# add_executable(testGreensFunction testGreensFunction.cxx)
# target_link_libraries(testGreensFunction eisvogel)
