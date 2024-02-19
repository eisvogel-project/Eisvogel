set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

add_executable(testDistributedNDArray testDistributedNDArray.cxx)
target_link_libraries(testDistributedNDArray eisvogel)

add_executable(testDistributedNDArrayDiagnostics testDistributedNDArrayDiagnostics.cxx)
target_link_libraries(testDistributedNDArrayDiagnostics eisvogel)

add_executable(testDistributedWeightingField testDistributedWeightingField.cxx)
target_link_libraries(testDistributedWeightingField eisvogel)

add_executable(testSerialization testSerialization.cxx)
target_link_libraries(testSerialization eisvogel)
