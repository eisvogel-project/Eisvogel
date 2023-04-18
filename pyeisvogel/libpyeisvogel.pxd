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

cdef extern from "Eisvogel/Common.hh":
    ctypedef float scalar_t
