#ifndef __NDARRAY_OPERATIONS_HH
#define __NDARRAY_OPERATIONS_HH

#include "DenseNDArray.hh"
#include "Eisvogel/IteratorUtils.hh"
#include <cmath>

namespace NDArrayOps {

  template <class T, std::size_t dims>
  DenseNDArray<T, dims> concatenate(const std::vector<DenseNDArray<T, dims>>& arrs, std::size_t axis) {

    if(axis >= dims) {
      [[unlikely]];
      throw std::runtime_error("Error: 'axis' out of bounds");
    }

    if(arrs.size() == 0) {
      [[unlikely]];
      throw std::runtime_error("Error: nothing to concatenate!");
    }

    if(arrs.size() == 1) {
      [[unlikely]];
      return arrs[0];
    }

    IndexVector output_shape = arrs[0].shape();

    // Check if shapes are compatible to allow concatenation: shapes of all arrays must be identical in all directions except `axis`
    for(std::size_t arr_ind = 1; arr_ind < arrs.size(); arr_ind++) {
      for(std::size_t dim = 0; dim < dims; dim++) {
	if(dim == axis) {
	  output_shape(dim) = output_shape(dim) + arrs[arr_ind].shape(dim);
	  continue;
	}
	if(output_shape(dim) != arrs[arr_ind].shape(dim)) {
	  throw std::runtime_error("Error: shapes incompatible; cannot concatenate");
	}
      }
    }

    std::cout << "HHH concatenating " << arrs.size() << " arrays" << std::endl;
    
    std::cout << "HHHH have concatenated output with dimension of" << std::endl;
    output_shape.print();
    std::cout << "HHHH" << std::endl;
    
    // --- this is just a crutch for now until we have fixed-size vectors
    std::array<std::size_t, dims> output_shape_crutch;
    std::copy(std::begin(output_shape), std::end(output_shape), std::begin(output_shape_crutch));
    // ---------
    
    DenseNDArray<T, dims> retval(output_shape_crutch, 0.0);

    // start assembling concatenated array
    std::size_t tmp_cnt = 0;
    IndexVector offset(dims, 0);
    for(auto& cur_arr : arrs) {

      // Migrate contents of `cur_arr` into output array
      IndexVector start_inds(dims, 0);
      IndexVector end_inds = cur_arr.shape();
      for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
	IndexVector cur_ind = cnt.index();
	IndexVector new_ind = cur_ind + offset;
	retval(new_ind) = cur_arr(cur_ind);
      }

      offset(axis) = offset(axis) + cur_arr.shape(axis);

      std::cout << "just migrated array #" << tmp_cnt << std::endl;
      tmp_cnt++;
    }

    return retval;
  }
  
  template <class T, std::size_t dims>
  DenseNDArray<T, dims> concatenate(const DenseNDArray<T, dims>& arr_1, const DenseNDArray<T, dims>& arr_2, std::size_t axis) {

    // std::cout << " ---> START CONCATENATE <---" << std::endl;
    
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

    // std::cout << " ---> MIGRATE ARR_1 <---" << std::endl;
    
    // Migrate contents of arr_1
    {
      IndexVector start_inds(dims, 0);
      IndexVector end_inds = shape_1;
      for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
	IndexVector cur_ind = cnt.index();
	retval(cur_ind) = arr_1(cur_ind);
      }
    }

    // std::cout << " ---> MIGRATE ARR_2 <---" << std::endl;
    
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

    // std::cout << " ---> FINISH CONCATENATE <---" << std::endl;
    
    return retval;
  }

  template <class T, std::size_t dims, template <class, std::size_t> class ArrayT>
  DenseNDArray<T, dims> range(const ArrayT<T, dims>& arr, const IndexVector& start_inds, const IndexVector& stop_inds) {

    if((start_inds.size() != dims) || (stop_inds.size() != dims)) {
      throw std::runtime_error("Error: not a possible range!");
    }

    IndexVector range_shape = stop_inds - start_inds;

    // --- this is just a crutch for now until we have fixed-size vectors
    std::array<std::size_t, dims> range_shape_crutch;
    std::copy(std::begin(range_shape), std::end(range_shape), std::begin(range_shape_crutch));
    // ---------

    DenseNDArray<T, dims> retval(range_shape_crutch, T());    
    for(IndexCounter cnt(start_inds, stop_inds); cnt.running(); ++cnt) {
      IndexVector cur_ind = cnt.index();
      IndexVector range_ind = cur_ind - start_inds;
      retval(range_ind) = arr(cur_ind);
    }
   
    return retval;
  }
  
  template<class T, std::size_t dims, template <class, std::size_t> class ArrayT>
  DenseNDArray<T, dims> min(const ArrayT<T, dims>& a, const ArrayT<T, dims>& b) {

    DenseNDArray<T, dims> retval(a.shape(), T());
      
    IndexVector start_inds(dims, 0);
    IndexVector end_inds = a.shape();
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
      IndexVector cur_ind = cnt.index();
      retval(cur_ind) = std::min(a(cur_ind), b(cur_ind));
    }    
    return retval;
  }

  template<class T, std::size_t dims>
  std::size_t number_nonzero_elems(const DenseNDArray<T, dims>& a) {
    std::size_t nonzero_cnt = 0;

    IndexVector start_inds(dims, 0);
    IndexVector end_inds = a.shape();
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
      IndexVector cur_ind = cnt.index();     
      if(a(cur_ind) != 0.0) {
	nonzero_cnt++;
      }
    }   
    return nonzero_cnt;
  }
  
}; // end namespace
  
#endif
