set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -funroll-loops -march=native")

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)

add_library(greenstor SHARED
            H5Utils.cxx)
target_include_directories(greenstor PUBLIC ${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
target_link_libraries(greenstor PUBLIC ${HDF5_LIBRARIES} MPI::MPI_CXX)
target_compile_options(greenstor PUBLIC ${HDF5_CFLAGS_OTHER})
set_target_properties(greenstor PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/greenstor)

add_library(eisvogel SHARED
            WeightingFieldUtils.cxx
            Kernels.cxx
            Interpolator.cxx
            Integrator.cxx
            Trajectory.cxx
            MathUtils.cxx
            SignalExport.cxx
            SignalCalculator.cxx)
target_include_directories(eisvogel PUBLIC ${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS})
set_target_properties(eisvogel PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel)

if(BUILD_WFC)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(MEEP REQUIRED meep)

    add_library(eisvogel_wfc SHARED
                WeightingFieldCalculator.cxx
                Geometry.cxx
                Antenna.cxx)
    target_link_directories(eisvogel_wfc PUBLIC ${MEEP_LIBRARY_DIRS} ${HDF5_LIBRARY_DIRS})
    target_link_libraries(eisvogel_wfc PUBLIC eisvogel greenstor ${MEEP_LIBRARIES} ${HDF5_LIBRARIES} MPI::MPI_CXX)
    target_include_directories(eisvogel_wfc PUBLIC ${PROJECT_SOURCE_DIR}/include ${MEEP_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
    target_compile_options(eisvogel_wfc PUBLIC ${MEEP_CFLAGS_OTHER} ${HDF5_CFLAGS_OTHER})

    set_target_properties(eisvogel_wfc PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/eisvogel_wfc)
endif(BUILD_WFC)
