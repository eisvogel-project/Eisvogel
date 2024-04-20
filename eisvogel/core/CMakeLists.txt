find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)
pkg_check_modules(UUID REQUIRED uuid)

add_library(eisvogel INTERFACE)
target_include_directories(eisvogel INTERFACE
	${PROJECT_SOURCE_DIR}/eisvogel/core
	${PROJECT_SOURCE_DIR}/eisvogel/core/impl)
target_link_libraries(eisvogel INTERFACE
	${UUID_LIBRARIES} ${TBB_LIBRARIES})
