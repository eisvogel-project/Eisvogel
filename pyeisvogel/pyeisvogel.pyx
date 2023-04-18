from pyeisvogel.libpyeisvogel cimport *
from libcpp.utility cimport move
import os

from pyeisvogel cimport ccoordutils
cdef class CoordVector:
     cdef ccoordutils.CoordVector* c_vec

     @staticmethod
     cdef __c_MakeCoordVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z):
         cdef vector[scalar_t] vec
         vec.push_back(t)
         vec.push_back(z)
         vec.push_back(x)
         vec.push_back(y)

         cdef CoordVector cv = CoordVector.__new__(CoordVector)
         cv.c_vec = new ccoordutils.CoordVector(move(vec))
         return cv

     @staticmethod
     def MakeCoordVectorTXYZ(t, x, y, z):
         return CoordVector.__c_MakeCoordVectorTXYZ(t, x, y, z)

from pyeisvogel cimport ccurrent
cdef class Current0D:
     cdef ccurrent.Current0D* c_current

     def __init__(self, points, charges):

         

         pass

from pyeisvogel cimport csignalcalculator
cdef class SignalCalculator:
    cdef csignalcalculator.SignalCalculator* c_calc

    def __init__(self, geometry_path):
        self.c_calc = new csignalcalculator.SignalCalculator(geometry_path.encode("utf-8"))

    def ComputeSignal(track, t_sig):
        return 0.0
