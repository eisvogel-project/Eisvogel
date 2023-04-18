from pyeisvogel.libpyeisvogel cimport *

cdef extern from "Eisvogel/CoordUtils.hh":
     cdef cppclass CoordVector:
          CoordVector(vector[scalar_t]&& data)
