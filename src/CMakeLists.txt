set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

find_package(PkgConfig REQUIRED)
pkg_check_modules(TBB REQUIRED tbb)
pkg_check_modules(UUID REQUIRED uuid)

add_library(eisvogel SHARED
	    IteratorUtils.cxx
	    WeightingFieldUtilsOld.cxx
            Trajectory.cxx
            MathUtils.cxx
            SignalExport.cxx
            SignalCalculatorOld.cxx
	    FieldStorage.cxx)
target_include_directories(eisvogel
	PUBLIC "${PROJECT_SOURCE_DIR}/include"
	PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(eisvogel PUBLIC ${UUID_LIBRARIES} ${TBB_LIBRARIES})
set_target_properties(eisvogel PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel)

if(BUILD_WFC)
    find_package(MPI REQUIRED)
    pkg_check_modules(MEEP REQUIRED meep)

    add_library(eisvogel_wfc_old SHARED
                CylindricalWeightingFieldCalculator.cxx
                GeometryOld.cxx
                AntennaOld.cxx)
    target_link_directories(eisvogel_wfc_old PUBLIC ${MEEP_LIBRARY_DIRS})
    target_link_libraries(eisvogel_wfc_old PUBLIC eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel_wfc_old PUBLIC ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS})
    target_compile_options(eisvogel_wfc_old PUBLIC ${MEEP_CFLAGS_OTHER})

    set_target_properties(eisvogel_wfc_old PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel_wfc_old)


    add_library(eisvogel_wfc SHARED
                Geometry.cxx)
    target_link_directories(eisvogel_wfc PUBLIC ${MEEP_LIBRARY_DIRS})
    target_link_libraries(eisvogel_wfc PUBLIC eisvogel ${MEEP_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel_wfc PUBLIC ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS})
    target_compile_options(eisvogel_wfc PUBLIC ${MEEP_CFLAGS_OTHER})

    set_target_properties(eisvogel_wfc PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel_wfc)

endif(BUILD_WFC)
