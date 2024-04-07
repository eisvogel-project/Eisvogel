namespace Interpolation::Kernel {

  void Linear::kern(scalar_t i_frac, scalar_t* out) {

    assert((i_frac >= 0) && (i_frac < 1));

    *(out+0) = 1 - i_frac;  // weight for sample at `i_int`
    *(out+1) = i_frac;      // weight for sample at `i_int + 1`
  }
  
  void Keys::kern(scalar_t i_frac, scalar_t* out) {

    // The total `scalar_t`-valued interpolation target `i` is split up into
    // `i = i_frac + i_int`, where `i_frac` is the fractional part and `i_int` the integer part
    
    // `i_frac` is supposed to be the fractional part for the interpolation
    assert((i_frac >= 0) && (i_frac < 1));

    // Evaluate the kernel at the four participating positions
    *(out+0) = i_frac * (-0.5 + (1.0 - 0.5 * i_frac) * i_frac);    // weight for sample at `i_int - 1`
    *(out+1) = 1.0 + i_frac * i_frac * (-2.5 + 1.5 * i_frac);      // weight for sample at `i_int`
    *(out+2) = i_frac * (0.5 + (2.0 - 1.5 * i_frac) * i_frac);     // weight for sample at `i_int + 1`
    *(out+3) = (-0.5 + 0.5 * i_frac) * i_frac * i_frac;            // weight for sample at `i_int + 2`
  }
}

namespace Interpolation {

  template <class KernelT, typename T, std::size_t dims, std::size_t vec_dims>
  void interpolate(const NDVecArray<T, dims, vec_dims>& arr, NDVecArray<T, 1, vec_dims>& retval,
		   const Vector<scalar_t, dims-1>& outer_inds,
		   scalar_t inner_ind_start, scalar_t inner_ind_end, scalar_t inner_ind_delta) {

    // Note: requirement will be loosened later
    static_assert(dims == 3);

    
    
    // compute int + fractional part of outer_inds
    
    // compute dims-1 kernel array for outer inds
    
    // iterate over interpolation range needed for outer_inds
    
    // inner loop over contiguous `inner_ind` direction
    
    // call 1d interpolations, weight with outer_ind kernel matrix and add to retval
    
  }
}
