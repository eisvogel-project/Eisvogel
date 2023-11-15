set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(02shower_simulation2d 02calculate_shower_emission2d.cpp)
target_link_libraries(02shower_simulation2d eisvogel)
target_link_libraries(02shower_simulation2d ShowerCreator)


target_include_directories(02shower_simulation2d PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_creator"
                          )
configure_file(${PROJECT_SOURCE_DIR}/examples/shower_profile/plot_results.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY
)