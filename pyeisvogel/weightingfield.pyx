from pyeisvogel.libpyeisvogel cimport *
from pyeisvogel cimport cweightingfield
from pyeisvogel.cserialization cimport Serializer
import os

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
