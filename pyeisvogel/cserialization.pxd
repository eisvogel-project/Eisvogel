cdef extern from "<fstream>" namespace "std":
     cdef cppclass fstream:
          pass

cdef extern from "Eisvogel/Serialization.hh" namespace "stor":
     cdef cppclass Serializer:
          Serializer(fstream& stream)
          T deserialize[T]()
