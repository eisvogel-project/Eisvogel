from python.libeisvogel cimport *
from python.ccoordutils cimport *
from libcpp.vector cimport vector

cdef extern from "Eisvogel/Current0D.hh":
     cdef cppclass Current0D:
          Current0D() except +
          Current0D(vector[CoordVector]&& points, vector[scalar_t]&& charges)
