set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(surface_sweep surface_sweep.cxx)
target_link_libraries(surface_sweep eisvogel)
