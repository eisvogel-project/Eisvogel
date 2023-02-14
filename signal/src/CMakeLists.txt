add_library(signal DenseWeightingFieldProvider.cxx)

target_include_directories(signal PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
