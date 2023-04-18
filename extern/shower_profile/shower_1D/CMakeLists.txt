add_library(Shower1D shower_1D.cpp)

target_include_directories(Shower1D PUBLIC
                          "${PROJECT_SOURCE_DIR}/environment"
                          "${PROJECT_SOURCE_DIR}/charge_excess_profile"
                          )