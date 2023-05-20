from python.libeisvogel cimport *
from python.ccurrent cimport Current0D, SparseCurrentDensity3D
from libcpp.string cimport string

cdef extern from "Eisvogel/SignalCalculator.hh":
     cdef cppclass SignalCalculator:
          SignalCalculator(const string& geometry_path)
          scalar_t ComputeSignal(const Current0D& track, scalar_t t_sig)
          scalar_t ComputeSignal(const SparseCurrentDensity3D& current_density, scalar_t t_sig)
