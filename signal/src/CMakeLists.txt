add_library(signal WeightingFieldUtils.cxx Kernels.cxx)

target_include_directories(signal PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
