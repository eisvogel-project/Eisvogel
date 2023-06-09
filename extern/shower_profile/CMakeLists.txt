cmake_minimum_required( VERSION 3.20 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_library(Shower1D SHARED shower_1D/shower_1D.cpp)
add_library(Shower2D SHARED shower_2D/shower_2D.cpp)
add_library(Environment SHARED environment/ice_profile.cpp)
add_library(ShowerCreator SHARED shower_creator/shower_creator.cpp)

target_link_libraries(ShowerCreator Shower1D)
target_link_libraries(ShowerCreator Shower2D)
target_link_libraries(ShowerCreator Environment)
target_link_libraries(Shower1D PUBLIC Environment)
target_link_libraries(Shower2D PUBLIC Environment)

target_include_directories(Shower1D PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/environment"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/charge_excess_profile"
                          "${PROJECT_SOURCE_DIR}/extern/utilities"
                          "${PROJECT_SOURCE_DIR}/extern/constants"
                          "${PROJECT_SOURCE_DIR}/include/"
                          )

target_include_directories(ShowerCreator PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_1D"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_2D"
                          "${PROJECT_SOURCE_DIR}/extern/utilities"
                          )
target_include_directories(Environment PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/utilities"
                          )

target_include_directories(Shower2D PUBLIC
                        "${PROJECT_SOURCE_DIR}/extern/shower_profile/environment"
                        "${PROJECT_SOURCE_DIR}/extern/shower_profile/charge_excess_profile"
                        "${PROJECT_SOURCE_DIR}/extern/utilities"
                        "${PROJECT_SOURCE_DIR}/extern/constants"
                        "${PROJECT_SOURCE_DIR}/include/"
                        )