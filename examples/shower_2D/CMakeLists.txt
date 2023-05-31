set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(shower_simulation_2D shower_profile_2D.cpp)
target_link_libraries(shower_simulation_2D eisvogel)
target_link_libraries(shower_simulation_2D ShowerCreator)

target_include_directories(shower_simulation_2D PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_creator"
                          )