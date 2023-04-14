cdef extern from "Common.hh":
    ctypedef float scalar_t

cdef extern from "Kernels.hh":
     cdef cppclass KeysCubicInterpolationKernel:
          KeysCubicInterpolationKernel() except +
          size_t Support() const
          scalar_t operator()(scalar_t arg) const