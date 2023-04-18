cdef extern from "WeightingField.hh":
     cdef cppclass WeightingField:
          WeightingField(WeightingField&& other)