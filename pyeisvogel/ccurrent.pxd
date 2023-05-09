from pyeisvogel.libpyeisvogel cimport *
from pyeisvogel.ccoordutils cimport *

cdef extern from "Eisvogel/Current0D.hh":
     cdef cppclass Current0D:
          Current0D() except +
          Current0D(vector[CoordVector]&& points, vector[scalar_t]&& charges)
