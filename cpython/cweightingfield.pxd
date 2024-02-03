from cpython.libeisvogel cimport *
from cpython.ccoordutils cimport CoordVector
from libcpp.string cimport string

cdef extern from "Eisvogel/WeightingField.hh":
     cdef cppclass CylindricalWeightingField:
          CylindricalWeightingField(string wf_path);
          scalar_t E_r(CoordVector& pos);
          scalar_t E_z(CoordVector& pos);
          scalar_t GetStartCoords(size_t dim) const;
          scalar_t GetEndCoords(size_t dim) const;  
