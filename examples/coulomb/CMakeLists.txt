set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(coulomb coulomb.cxx)
target_link_libraries(coulomb eisvogel)
