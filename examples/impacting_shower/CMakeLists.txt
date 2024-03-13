set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(impacting_shower impacting_shower.cxx)
target_link_libraries(impacting_shower eisvogel)
