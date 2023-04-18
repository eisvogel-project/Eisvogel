cimport ckernels

cdef class KeysCubicInterpolationKernel:
     cdef ckernels.KeysCubicInterpolationKernel c_kernel

     def evaluate(self, x):
         return self.c_kernel(x)