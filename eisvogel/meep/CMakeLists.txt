find_package(MPI REQUIRED)
pkg_check_modules(MEEP REQUIRED meep)

add_library(eisvogel-meep INTERFACE)
target_include_directories(eisvogel-meep INTERFACE
       ${PROJECT_SOURCE_DIR}/eisvogel/meep
       ${PROJECT_SOURCE_DIR}/eisvogel/meep/impl
       ${MEEP_INCLUDE_DIRS})
target_link_directories(eisvogel-meep INTERFACE ${MEEP_LIBRARY_DIRS})
target_link_libraries(eisvogel-meep INTERFACE eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
target_compile_options(eisvogel-meep INTERFACE ${MEEP_CFLAGS_OTHER})
