cdef extern from "Eisvogel/WeightingField.hh":
     cdef cppclass WeightingField:
          WeightingField(WeightingField&& other)