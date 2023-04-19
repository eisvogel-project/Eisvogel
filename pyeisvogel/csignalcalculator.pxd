from pyeisvogel.libpyeisvogel cimport *
from pyeisvogel.ccurrent cimport Current0D

cdef extern from "Eisvogel/SignalCalculator.hh":
     cdef cppclass SignalCalculator:
          SignalCalculator(const string& geometry_path)
          scalar_t ComputeSignal(const Current0D& track, scalar_t t_sig)
