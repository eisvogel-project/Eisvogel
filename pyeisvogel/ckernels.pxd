from pyeisvogel.libpyeisvogel cimport *

cdef extern from "Eisvogel/Kernels.hh":
     cdef cppclass KeysCubicInterpolationKernel:
          KeysCubicInterpolationKernel() except +
          size_t Support() const
          scalar_t operator()(scalar_t arg) const
