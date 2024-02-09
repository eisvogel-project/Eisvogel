#ifndef __NDARRAY_OPERATIONS_HH
#define __NDARRAY_OPERATIONS_HH

#include "DenseNDArray.hh"
#include "Eisvogel/IteratorUtils.hh"

namespace NDArrayOps {

  template <class T, std::size_t dims>
  DenseNDArray<T, dims> concatenate(const DenseNDArray<T, dims>& arr_1, const DenseNDArray<T, dims>& arr_2, std::size_t axis) {

    if(axis >= dims) {
      throw std::runtime_error("Error: 'axis' out of bounds");
    }
    
    IndexVector shape_1 = arr_1.shape();
    IndexVector shape_2 = arr_2.shape();
    
    // Check if shapes are compatible to allow concatenation: shapes must be identical in all directions except `axis`
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(ind == axis) {
	continue;
      }
      if(shape_1(ind) != shape_2(ind)) {
	throw std::runtime_error("Error: shapes incompatible; cannot concatenate");
      }
    }
    
    IndexVector final_shape = shape_1;
    final_shape(axis) = shape_1(axis) + shape_2(axis);

    // --- this is just a crutch for now until we have fixed-size vectors
    std::array<std::size_t, dims> final_shape_crutch;
    std::copy(std::begin(final_shape), std::end(final_shape), std::begin(final_shape_crutch));
    // ---------
    
    DenseNDArray<T, dims> retval(final_shape_crutch, 0.0);

    // Migrate contents of arr_1
    {
      IndexVector start_inds(dims, 0);
      IndexVector end_inds = shape_1;
      for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
	IndexVector cur_ind = cnt.index();
	retval(cur_ind) = arr_1(cur_ind);
      }
    }

    // Migrate contents of arr_2
    {
      // `arr_2` is offset only along the direction the concatenation is performed
      IndexVector offset(dims, 0);
      offset(axis) = shape_1(axis);
      
      IndexVector start_inds(dims, 0);
      IndexVector end_inds = shape_2;
      for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
	IndexVector cur_ind = cnt.index();
	IndexVector new_ind = cur_ind + offset;
	retval(new_ind) = arr_2(cur_ind);
      }
    }
    
    return retval;
  }

  template <class T, std::size_t dims>
  DenseNDArray<T, dims> range(const DenseNDArray<T, dims>& arr, const IndexVector& start_inds, const IndexVector& stop_inds) {

    if((start_inds.size() != dims) || (stop_inds.size() != dims)) {
      throw std::runtime_error("Error: not a possible range!");
    }

    IndexVector range_shape = stop_inds - start_inds;

    // --- this is just a crutch for now until we have fixed-size vectors
    std::array<std::size_t, dims> range_shape_crutch;
    std::copy(std::begin(range_shape), std::end(range_shape), std::begin(range_shape_crutch));
    // ---------

    DenseNDArray<float, dims> retval(range_shape_crutch, 0.0);    
    for(IndexCounter cnt(start_inds, stop_inds); cnt.running(); ++cnt) {
      IndexVector cur_ind = cnt.index();
      IndexVector range_ind = cur_ind - start_inds;
      retval(range_ind) = arr(cur_ind);
    }
   
    return retval;
  }  
  
}; // end namespace
  
#endif
