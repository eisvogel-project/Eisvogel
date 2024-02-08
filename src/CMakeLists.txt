set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

find_package(PkgConfig REQUIRED)
pkg_check_modules(UUID REQUIRED uuid)

add_library(eisvogel SHARED
            WeightingFieldUtils.cxx
	    IteratorUtils.cxx
            Trajectory.cxx
            MathUtils.cxx
            SignalExport.cxx
            SignalCalculator.cxx
	    FieldStorage.cxx)
target_include_directories(eisvogel
	PUBLIC "${PROJECT_SOURCE_DIR}/include" ${MPI_INCLUDE_PATH}
	PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(eisvogel PUBLIC MPI::MPI_CXX ${UUID_LIBRARIES})
set_target_properties(eisvogel PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel)

if(BUILD_WFC)
    find_package(MPI REQUIRED)
    pkg_check_modules(MEEP REQUIRED meep)

    add_library(eisvogel_wfc SHARED
                CylindricalWeightingFieldCalculator.cxx
		MPIUtils.cxx
                Geometry.cxx
                Antenna.cxx)
    target_link_directories(eisvogel_wfc PUBLIC ${MEEP_LIBRARY_DIRS})
    target_link_libraries(eisvogel_wfc PUBLIC eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel_wfc PUBLIC ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS})
    target_compile_options(eisvogel_wfc PUBLIC ${MEEP_CFLAGS_OTHER})

    set_target_properties(eisvogel_wfc PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel_wfc)
endif(BUILD_WFC)
