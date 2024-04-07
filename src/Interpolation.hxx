namespace Interpolation {

  template <class KernelT, typename T, std::size_t dims, std::size_t vec_dims>
  void interpolate(const NDVecArray<T, dims, vec_dims>& arr, NDVecArray<T, 1, vec_dims>& retval,
		   const Vector<scalar_t, dims-1>& outer_inds,
		   scalar_t inner_ind_start, scalar_t inner_ind_end, scalar_t inner_ind_delta) {
    
    // Note: requirement on `dims == 3` will be loosened later

    std::cout << KernelT::support << std::endl;
    
    std::cout << arr.m_data -> size() << std::endl;
    
    // compute int + fractional part of outer_inds
    
    // compute dims-1 kernel array for outer inds
    
    // iterate over interpolation range needed for outer_inds
    
    // inner loop over contiguous `inner_ind` direction
    
    // call 1d interpolations, weight with outer_ind kernel matrix and add to retval
    
  }
}
