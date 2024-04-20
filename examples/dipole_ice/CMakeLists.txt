cmake_minimum_required( VERSION 3.20 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(dipole_ice dipole_ice.cxx)
target_link_libraries(dipole_ice eisvogel_gfc)
