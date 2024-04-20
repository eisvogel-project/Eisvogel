set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)
pkg_check_modules(UUID REQUIRED uuid)

add_library(eisvogel SHARED
            SignalExport.cxx
            SignalCalculator.cxx)
target_include_directories(eisvogel
	PUBLIC "${PROJECT_SOURCE_DIR}/include"
	PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(eisvogel PUBLIC ${UUID_LIBRARIES} ${TBB_LIBRARIES})
set_target_properties(eisvogel PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel)

if(BUILD_GFC)

    find_package(MPI REQUIRED)
    pkg_check_modules(MEEP REQUIRED meep)

    add_library(eisvogel_gfc SHARED
                Geometry.cxx
		Antenna.cxx
		CylindricalGreensFunctionCalculator.cxx)
    target_link_directories(eisvogel_gfc PUBLIC ${MEEP_LIBRARY_DIRS})
    target_link_libraries(eisvogel_gfc PUBLIC eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel_gfc PUBLIC ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS})
    target_compile_options(eisvogel_gfc PUBLIC ${MEEP_CFLAGS_OTHER})

    set_target_properties(eisvogel_gfc PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel_gfc)

endif(BUILD_GFC)
