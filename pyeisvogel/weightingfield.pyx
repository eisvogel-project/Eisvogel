from pyeisvogel cimport cweightingfield
from pyeisvogel.cserialization cimport Serializer
import os

cdef extern from "<iostream>":
     cdef cppclass openmode:
          pass
     cdef openmode binary_in "std::ios_base::binary | std::ios_base::in"

cdef extern from "<iostream>" namespace "std":
     cdef cppclass iostream:
          pass

cdef extern from "<fstream>" namespace "std":
     cdef cppclass fstream(iostream):
          void open(const char*, openmode)

cdef extern from "<string>" namespace "std":
     cdef cppclass string:
          char* c_str()
          string(char*)

cdef class WeightingField:
     cdef cweightingfield.WeightingField* c_wf

     def __cinit__(self):
         pass

     @staticmethod
     cdef WeightingField from_path_c(string& path):
         cdef fstream ifs
         ifs.open(path.c_str(), binary_in)

         cdef Serializer* ser = new Serializer(ifs)
         cdef WeightingField wf = WeightingField.__new__(WeightingField)
         wf.c_wf = new cweightingfield.WeightingField(ser.deserialize[cweightingfield.WeightingField]())
         return wf

     @staticmethod
     def from_path(path):
         if not os.path.exists(path):
             raise RuntimeError("Error: file does not exist!")

         return WeightingField.from_path_c(path.encode('utf-8'))
