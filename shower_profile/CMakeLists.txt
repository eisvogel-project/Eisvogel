cmake_minimum_required(VERSION 3.10)

project(Shower)
add_subdirectory(environment)
add_subdirectory(charge_excess_profile)
add_subdirectory(shower_1D)
add_subdirectory(shower_creator)
add_executable(Shower main.cpp)

target_link_libraries(Shower Shower1D)
target_link_libraries(Shower ShowerCreator)
target_link_libraries(Shower ROOT::RIO ROOT::Net)
target_include_directories(Shower PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/shower_1D"
                          "${PROJECT_SOURCE_DIR}/shower_creator"
                          )
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Shower ${ROOT_LIBRARIES})