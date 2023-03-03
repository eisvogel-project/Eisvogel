add_library(signal WeightingFieldUtils.cxx Kernels.cxx Interpolator.cxx Integrator.cxx Trajectory.cxx MathUtils.cxx SignalExport.cxx)

target_include_directories(signal PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
