add_library(ShowerCreator shower_creator.cpp)

target_link_libraries(ShowerCreator Shower1D)
target_link_libraries(ShowerCreator Shower2D)
target_link_libraries(ShowerCreator Environment)

target_include_directories(ShowerCreator PUBLIC
                          "${PROJECT_SOURCE_DIR}/environment"
                          "${PROJECT_SOURCE_DIR}/charge_excess_profile"
                          "${PROJECT_SOURCE_DIR}/shower_1D"
                          "${PROJECT_SOURCE_DIR}/shower_2D"
                          )