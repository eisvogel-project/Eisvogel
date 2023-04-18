cimport cweightingfield

cdef class WeightingField:
     cdef cweightingfield.WeightingField* c_wf

     def __cinit__(self):
         pass

     # @staticmethod
     # cdef WeightingField from_path(cweightingfield.WeightingField* path):
     #     cdef WeightingField wf = WeightingField.__new__(WeightingField)
     #     wf.c_wf = path
     #     return wf