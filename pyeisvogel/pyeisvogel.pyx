from pyeisvogel.libpyeisvogel cimport *
from libcpp.utility cimport move
from libcpp.memory cimport unique_ptr, make_unique
import os

from pyeisvogel cimport ccoordutils
cdef class CoordVector:
    cdef unique_ptr[ccoordutils.CoordVector] c_vec

    @staticmethod
    cdef __c_FromTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z):
        cdef CoordVector vec = CoordVector.__new__(CoordVector)
        vec.c_vec = make_unique[ccoordutils.CoordVector](ccoordutils.MakeCoordVectorTXYZ(t, x, y, z))
        return vec

    @staticmethod
    def FromTXYZ(t, x, y, z):
        return CoordVector.__c_FromTXYZ(t, x, y, z)

from pyeisvogel cimport ccurrent
cdef class Current0D:
    cdef ccurrent.Current0D* c_current
    
from pyeisvogel cimport csignalcalculator
cdef class SignalCalculator:
    cdef csignalcalculator.SignalCalculator* c_calc

    def __init__(self, geometry_path):
        self.c_calc = new csignalcalculator.SignalCalculator(geometry_path.encode("utf-8"))

    def ComputeSignal(track, t_sig):
        return 0.0
