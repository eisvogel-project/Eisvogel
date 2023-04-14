add_library(signal SHARED WeightingFieldUtils.cxx Kernels.cxx Interpolator.cxx Integrator.cxx Trajectory.cxx MathUtils.cxx SignalExport.cxx)

target_include_directories(signal PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))

set_target_properties(signal PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/signal)
