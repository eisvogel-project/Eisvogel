from cython.operator import dereference
from python.libeisvogel cimport *
from libcpp.utility cimport move
from libcpp.memory cimport unique_ptr, make_unique
from libcpp.vector cimport vector
from libcpp.string cimport string
import os

from python cimport ccoordutils
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

from python cimport ccurrent
cdef class Current0D:
    cdef unique_ptr[ccurrent.Current0D] c_current

    @staticmethod
    def FromSegments(points, charges):
        cdef vector[ccoordutils.CoordVector] vec_points
        cdef vector[scalar_t] vec_charges
        cdef scalar_t charge
        cdef CoordVector point

        for point in points:
            vec_points.push_back(dereference(point.c_vec))

        for charge in charges:
            vec_charges.push_back(charge)

        cdef Current0D cur = Current0D.__new__(Current0D)
        cur.c_current = make_unique[ccurrent.Current0D](move(vec_points), move(vec_charges))
        return cur
    
from python cimport csignalcalculator
cdef class SignalCalculator:
    cdef csignalcalculator.SignalCalculator* c_calc

    def __init__(self, geometry_path):
        self.c_calc = new csignalcalculator.SignalCalculator(geometry_path.encode("utf-8"))

    def ComputeSignal(self, Current0D track, t_sig):
        cdef scalar_t signal
        signal = self.c_calc.ComputeSignal(dereference(track.c_current), t_sig)
        return signal

from python cimport cweightingfieldutils
cpdef CreateElectricDipoleWeightingField(string wf_path, CoordVector start_coords, CoordVector end_coords, 
                                         scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor):
    cweightingfieldutils.CreateElectricDipoleWeightingField(wf_path, 
                                                            dereference(start_coords.c_vec), dereference(end_coords.c_vec), 
                                                            tp, N, r_min, os_factor)
