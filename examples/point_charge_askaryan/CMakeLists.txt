set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(point_charge_askaryan point_charge_askaryan.cxx)
target_link_libraries(point_charge_askaryan eisvogel)
