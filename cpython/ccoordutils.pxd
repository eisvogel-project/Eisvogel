from cpython.libeisvogel cimport *
from libcpp.vector cimport vector

cdef extern from "Eisvogel/CoordUtils.hh":
     cdef cppclass CoordVector:
          CoordVector(vector[scalar_t]& data)
     cdef cppclass FieldVector:
          FieldVector(vector[scalar_t]& data)
     cdef cppclass DeltaVector:
          DeltaVector(vector[scalar_t]& data)
     cdef cppclass IndexVector:
          IndexVector(vector[size_t]& data)

cdef extern from "Eisvogel/CoordUtils.hh" namespace "CoordUtils":
     cdef CoordVector MakeCoordVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z)
     cdef CoordVector MakeCoordVectorTRZ(scalar_t t, scalar_t r, scalar_t z)

     cdef FieldVector MakeFieldVectorXYZ(scalar_t x, scalar_t y, scalar_t z)
     cdef DeltaVector MakeDeltaVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z)

     cdef IndexVector MakeIndexVectorTRZ(size_t ind_t, size_t ind_r, size_t ind_z)