add_library(Shower1D SHARED shower_1D/shower_1D.cpp)
target_include_directories(Shower1D PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/environment"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/charge_excess_profile"
                          )
