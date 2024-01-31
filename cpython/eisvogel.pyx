from cython.operator import dereference
from cpython.libeisvogel cimport *
from libcpp.utility cimport move
from libcpp.memory cimport unique_ptr, make_unique
from libcpp.vector cimport vector
from libcpp.string cimport string
import os

from cpython cimport ccoordutils
cdef class CoordVector:
    cdef unique_ptr[ccoordutils.CoordVector] c_vec

    @staticmethod
    cdef __c_FromTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z):
        cdef CoordVector vec = CoordVector.__new__(CoordVector)
        vec.c_vec = make_unique[ccoordutils.CoordVector](ccoordutils.MakeCoordVectorTXYZ(t, x, y, z))
        return vec

    @staticmethod
    cdef __c_FromTRZ(scalar_t t, scalar_t r, scalar_t z):
        cdef CoordVector vec = CoordVector.__new__(CoordVector)
        vec.c_vec = make_unique[ccoordutils.CoordVector](ccoordutils.MakeCoordVectorTRZ(t, r, z))
        return vec

    @staticmethod
    def FromTXYZ(t, x, y, z):
        return CoordVector.__c_FromTXYZ(t, x, y, z)

    @staticmethod
    def FromTRZ(t, r, z):
        return CoordVector.__c_FromTRZ(t, r, z)

cdef class FieldVector:
    cdef unique_ptr[ccoordutils.FieldVector] c_vec

    @staticmethod
    cdef __c_FromXYZ(scalar_t x, scalar_t y, scalar_t z):
        cdef FieldVector vec = FieldVector.__new__(FieldVector)
        vec.c_vec = make_unique[ccoordutils.FieldVector](ccoordutils.MakeFieldVectorXYZ(x, y, z))
        return vec

    @staticmethod
    def FromXYZ(x, y, z):
        return FieldVector.__c_FromXYZ(x, y, z)

cdef class DeltaVector:
    cdef unique_ptr[ccoordutils.DeltaVector] c_vec

    @staticmethod
    cdef __c_FromDeltaTXYZ(scalar_t delta_t, scalar_t delta_x, scalar_t delta_y, scalar_t delta_z):
        cdef DeltaVector vec = DeltaVector.__new__(DeltaVector)
        vec.c_vec = make_unique[ccoordutils.DeltaVector](ccoordutils.MakeDeltaVectorTXYZ(delta_t, delta_x, delta_y, delta_z))
        return vec

    @staticmethod
    def FromDeltaTXYZ(delta_t, delta_x, delta_y, delta_z):
        return DeltaVector.__c_FromDeltaTXYZ(delta_t, delta_x, delta_y, delta_z)

cdef class IndexVector:
    cdef unique_ptr[ccoordutils.IndexVector] c_ind

    @staticmethod
    def __c_FromTRZ(size_t ind_t, size_t ind_r, size_t ind_z):
        cdef IndexVector ind = IndexVector.__new__(IndexVector)
        ind.c_ind = make_unique[ccoordutils.IndexVector](ccoordutils.MakeIndexVectorTRZ(ind_t, ind_r, ind_z))
        return ind

    @staticmethod
    def FromTRZ(ind_t, ind_r, ind_z):
        return IndexVector.__c_FromTRZ(ind_t, ind_r, ind_z)

from cpython cimport ccurrent
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

cdef class SparseCurrentDensity3D:
    cdef ccurrent.SparseCurrentDensity3D* c_current_density

    def __init__(self, DeltaVector voxel_size):
        self.c_current_density = new ccurrent.SparseCurrentDensity3D(dereference(voxel_size.c_vec))

    def addCurrentElement(self, CoordVector point, FieldVector current_density):
        self.c_current_density.addCurrentElement(dereference(point.c_vec), dereference(current_density.c_vec))
    
from cpython cimport csignalcalculator
cdef class SignalCalculator:
    cdef csignalcalculator.SignalCalculator* c_calc

    def __init__(self, geometry_path):
        self.c_calc = new csignalcalculator.SignalCalculator(geometry_path.encode("utf-8"))

    def ComputeSignal(self, source, scalar_t t_sig):
        if isinstance(source, Current0D):
            return self.c_calc.ComputeSignal(dereference(
                (<Current0D?>(source)).c_current
            ), t_sig)
        elif isinstance(source, SparseCurrentDensity3D):
            return self.c_calc.ComputeSignal(dereference(
                (<SparseCurrentDensity3D?>(source)).c_current_density
            ), t_sig)
        else:
            raise RuntimeError("Unknown source type")

from cpython cimport cdistributedweightingfield
cdef class DistributedWeightingField:
    cdef cdistributedweightingfield.DistributedWeightingField* c_dwf

    def __init__(self, wf_path):
        self.c_dwf = new cdistributedweightingfield.DistributedWeightingField(wf_path.encode("utf-8"))

    def E_r(self, ind_trz):
        ind_t, ind_r, ind_z = ind_trz
        cdef IndexVector ind = IndexVector.FromTRZ(ind_t, ind_r, ind_z)
        return self.c_dwf.E_r(dereference(ind.c_ind))

    def E_z(self, ind_trz):
        ind_t, ind_r, ind_z = ind_trz
        cdef IndexVector ind = IndexVector.FromTRZ(ind_t, ind_r, ind_z)
        return self.c_dwf.E_z(dereference(ind.c_ind))

    def E_phi(self, ind_trz):
        ind_t, ind_r, ind_z = ind_trz
        cdef IndexVector ind = IndexVector.FromTRZ(ind_t, ind_r, ind_z)
        return self.c_dwf.E_phi(dereference(ind.c_ind))

    def E_rzphi(self, ind_trz):
        ind_t, ind_r, ind_z = ind_trz
        cdef IndexVector ind = IndexVector.FromTRZ(ind_t, ind_r, ind_z)
        return [self.c_dwf.E_r(dereference(ind.c_ind)),  self.c_dwf.E_z(dereference(ind.c_ind)), self.c_dwf.E_phi(dereference(ind.c_ind))]

    def shape(self):
        return [self.c_dwf.shape(0), self.c_dwf.shape(2), self.c_dwf.shape(1)]

    def startInd(self):
        return [self.c_dwf.startInd(0), self.c_dwf.startInd(2), self.c_dwf.startInd(1)]

from cpython cimport cweightingfieldutils
cpdef CreateElectricDipoleWeightingField(wf_path, CoordVector start_coords, CoordVector end_coords, 
                                         scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor):
    cweightingfieldutils.CreateElectricDipoleWeightingField(wf_path.encode("utf-8"), 
                                                            dereference(start_coords.c_vec), 
                                                            dereference(end_coords.c_vec), 
                                                            tp, N, r_min, os_factor)
    
