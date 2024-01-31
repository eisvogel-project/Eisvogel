from cpython.libeisvogel cimport *
from cpython.ccoordutils cimport IndexVector
from libcpp.string cimport string

cdef extern from "Eisvogel/DistributedWeightingField.hh":
     cdef cppclass DistributedWeightingField:
          DistributedWeightingField(string wf_path);
          scalar_t E_r(IndexVector& ind);
          scalar_t E_z(IndexVector& ind);
          scalar_t E_phi(IndexVector& ind);
          size_t shape(size_t dim);
          size_t startInd(size_t dim);