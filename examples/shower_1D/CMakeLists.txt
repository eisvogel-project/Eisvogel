set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(shower_simulation_1D shower_profile_1D.cpp)
target_link_libraries(shower_simulation_1D eisvogel)
target_link_libraries(shower_simulation_1D ShowerCreator)

target_include_directories(shower_simulation_1D PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_creator"
                          )