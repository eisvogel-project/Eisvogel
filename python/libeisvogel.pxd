cdef extern from "<vector>" namespace "std":
     cdef cppclass vector[T]:
          vector() except +
          void push_back(T&) except +

cdef extern from "<string>" namespace "std":
     cdef cppclass string:
          char* c_str()
          string(char*)

cdef extern from "Eisvogel/Common.hh":
     ctypedef float scalar_t
