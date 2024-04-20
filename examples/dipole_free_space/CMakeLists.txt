set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(dipole_free_space dipole_free_space.cxx)
target_link_libraries(dipole_free_space eisvogel)
