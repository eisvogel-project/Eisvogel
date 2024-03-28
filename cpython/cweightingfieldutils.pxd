from cpython.libeisvogel cimport *
from cpython.ccoordutils cimport *
from libcpp.string cimport string

cdef extern from "Eisvogel/WeightingFieldUtilsOld.hh" namespace "WeightingFieldUtilsOld":
     void CreateElectricDipoleWeightingField(string wf_path, const CoordVector& start_coords, const CoordVector& end_coords,
                                             scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor)
