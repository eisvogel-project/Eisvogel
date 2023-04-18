cdef extern from "<iostream>" namespace "std":
     cdef cppclass iostream:
          pass

cdef extern from "Serialization.hh":
     cdef cppclass Serializer:
          Serializer(iostream& stream)
          T deserialize[T]()