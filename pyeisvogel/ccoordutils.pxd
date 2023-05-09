from pyeisvogel.libpyeisvogel cimport *

cdef extern from "Eisvogel/CoordUtils.hh":
     cdef cppclass CoordVector:
          CoordVector(vector[scalar_t]&& data)

cdef extern from "Eisvogel/CoordUtils.hh" namespace "CoordUtils":
     cdef CoordVector MakeCoordVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z)
