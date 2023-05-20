from python.libeisvogel cimport *
from python.ccoordutils cimport *
from libcpp.vector cimport vector

cdef extern from "Eisvogel/Current0D.hh":
     cdef cppclass Current0D:
          Current0D() except +
          Current0D(vector[CoordVector]&& points, vector[scalar_t]&& charges)

cdef extern from "Eisvogel/SparseCurrentDensity3D.hh":
     cdef cppclass SparseCurrentDensity3D:
          SparseCurrentDensity3D(const DeltaVector& voxel_size)
          void addCurrentElement(const CoordVector& point, const FieldVector& current_density)
