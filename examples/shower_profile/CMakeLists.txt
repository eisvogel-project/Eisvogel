set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(02shower_simulation 02calculate_shower_emission.cpp)
target_link_libraries(02shower_simulation eisvogel)
target_link_libraries(02shower_simulation ShowerCreator)

add_executable(01weighting_field 01generate_weighting_field.cpp)
target_link_libraries(01weighting_field eisvogel)
target_link_libraries(01weighting_field ShowerCreator)

add_executable(03_simple_shower 03_simple_shower.cpp)
target_link_libraries(03_simple_shower eisvogel)

target_include_directories(02shower_simulation PUBLIC
                          "${PROJECT_SOURCE_DIR}/extern/shower_profile/shower_creator"
                          )
configure_file(${PROJECT_SOURCE_DIR}/examples/shower_profile/plot_results.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY
)
