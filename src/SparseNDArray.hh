#ifndef __SPARSE_NDARRAY_HH
#define __SPARSE_NDARRAY_HH

#include "NDArray.hh"

template <class T, std::size_t dims>
class SparseNDArray : public NDArray<T, dims> {

public:
  using shape_t = typename NDArray<T, dims>::shape_t;
  
};

// Some type shortcuts
template <class T>
using SparseScalarField3D = SparseNDArray<T, 3>;

#endif
