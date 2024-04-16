set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(export_greens_function export_greens_function.cxx)
target_include_directories(export_greens_function
	PUBLIC "${PROJECT_SOURCE_DIR}/include"
	PUBLIC "${PROJECT_SOURCE_DIR}/src")
