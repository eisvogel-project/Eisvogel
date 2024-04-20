find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)
pkg_check_modules(UUID REQUIRED uuid)

add_library(eisvogel INTERFACE)
target_include_directories(eisvogel INTERFACE
	"${PROJECT_SOURCE_DIR}/include"
	"${PROJECT_SOURCE_DIR}/src")
target_link_libraries(eisvogel INTERFACE
	${UUID_LIBRARIES} ${TBB_LIBRARIES})

if(BUILD_MEEP)

    find_package(MPI REQUIRED)
    pkg_check_modules(MEEP REQUIRED meep)

    add_library(eisvogel-meep INTERFACE)
    target_link_directories(eisvogel-meep INTERFACE ${MEEP_LIBRARY_DIRS})
    target_link_libraries(eisvogel-meep INTERFACE eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel-meep INTERFACE ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS})
    target_compile_options(eisvogel-meep INTERFACE ${MEEP_CFLAGS_OTHER})

endif(BUILD_MEEP)
