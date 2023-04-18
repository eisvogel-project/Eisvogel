set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -msse2 -ftree-vectorizer-verbose=2 -funroll-loops -march=native -g -fno-omit-frame-pointer -fno-inline -fno-default-inline -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -msse2 -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

add_library(eisvogel SHARED WeightingFieldUtils.cxx Kernels.cxx Interpolator.cxx Integrator.cxx Trajectory.cxx MathUtils.cxx SignalExport.cxx SignalCalculator.cxx)

target_include_directories(eisvogel PUBLIC ${PROJECT_SOURCE_DIR}/include)

set_target_properties(eisvogel PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel)
