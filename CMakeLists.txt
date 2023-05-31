cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project(eisvogel)

add_subdirectory(cmake)
add_subdirectory(src)
add_subdirectory(cpython)
add_subdirectory(examples)
add_subdirectory(extern/shower_profile)

target_link_libraries(eisvogel PUBLIC Shower1D)
target_include_directories(eisvogel PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_1D"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_2D"
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/environment"
                          )