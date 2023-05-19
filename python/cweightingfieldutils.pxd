from python.libeisvogel cimport *
from python.ccordutils cimport *
from libcpp.string cimport string

cdef extern from "Eisvogel/WeightingFieldUtils.hh" namespace "WeightingFieldUtils":
     void CreateElectricDipoleWeightingField(string wf_path, const CoordVector& start_coords, const CoordVector& end_coords,
                                             scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor)
