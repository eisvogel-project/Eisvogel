#include <cmath>

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

  template <class KernelT, std::size_t vec_dims>
  void interpolate_1d(scalar_t* kernel, scalar_t* data, scalar_t* out) {
    
    constexpr std::size_t ksize = 2 * KernelT::support;

    scalar_t accum[vec_dims] = {0};
    for(std::size_t i = 0; i < ksize; i++) {
      for(std::size_t j = 0; j < vec_dims; j++) {       
	accum[j] += kernel[i] * data[vec_dims * i + j];
      }
    }

    for(std::size_t i = 0; i < vec_dims; i++) {
      out[i] = accum[i];
    }
  }
  
  template <class KernelT, typename T, std::size_t dims, std::size_t vec_dims>
  void interpolate(const NDVecArray<T, dims, vec_dims>& arr, NDVecArray<T, 1, vec_dims>& retval,
		   const Vector<scalar_t, dims-1>& outer_inds,
		   scalar_t inner_ind_start, scalar_t inner_ind_end, scalar_t inner_ind_delta) {

    // Note: requirement will be loosened later
    static_assert(dims == 3);

    constexpr std::size_t ksize = 2 * KernelT::support;

    // compute integer and fractional parts of outer indices
    std::size_t outer_int[dims - 1];
    scalar_t outer_frac[dims - 1];
    for(std::size_t i = 0; i < dims - 1; i++) {
      scalar_t cur_int = 0.0;
      outer_frac[i] = std::modf(outer_inds[i], &cur_int);
      outer_int[i] = (std::size_t)cur_int;
    }

    // evaluate interpolation kernel for the outer indices
    scalar_t outer_kernels[dims - 1][ksize];
    for(std::size_t i = 0; i < dims - 1; i++) {
      KernelT::kern(outer_frac[i], outer_kernels[i]);
    }

    // std::cout << "-- outer kernels --" << std::endl;
    // for(std::size_t i = 0; i < dims - 1; i++) {
    //   for(std::size_t j = 0; j < ksize; j++) {
    // 	std::cout << outer_kernels[i][j] << ", ";
    //   }
    //   std::cout << std::endl;
    // }    
    
    // compute direct product of outer kernel evaluations
    // TODO: this needs to be generalized to arbitrary dimensions
    scalar_t outer_weights[ksize][ksize];
    for(std::size_t i1 = 0; i1 < ksize; i1++) {
      for(std::size_t i2 = 0; i2 < ksize; i2++) {
	outer_weights[i1][i2] = outer_kernels[0][i1] * outer_kernels[1][i2];
      }
    }
    
    // iterate over interpolation range needed for outer_inds
    // TODO: this needs to be generalized to arbitrary dimensions

    // TODO: add outer batch loop here
    
    for(std::size_t i1 = outer_int[0] - 1; i1 < outer_int[0] - 1 + ksize; i1++) {
      for(std::size_t i2 = outer_int[1] - 1; i2 < outer_int[1] - 1 + ksize; i2++) {

	Vector<std::size_t, dims> column_ind(0);
	column_ind[0] = i1;
	column_ind[1] = i2;
	scalar_t* column_start = &(*(arr.m_data -> begin())) + arr.ComputeFlatInd(column_ind);

	std::size_t inner_elem = 0;
	for(scalar_t cur_inner_ind = inner_ind_start; cur_inner_ind < inner_ind_end; cur_inner_ind += inner_ind_delta) {

	  // std::cout << "on inner_elem = " << inner_elem << std::endl;
	  
	  scalar_t tmp = 0.0;
	  scalar_t inner_frac = std::modf(cur_inner_ind, &tmp);
	  std::size_t inner_int = (std::size_t)tmp;

	  // std::cout << "inner_frac = " << inner_frac << std::endl;
	  
	  scalar_t inner_interp[vec_dims] = {0};
	  scalar_t inner_kernel[ksize] = {0};
	  KernelT::kern(inner_frac, inner_kernel);

	  // std::cout << "inner_kernel = ";
	  // for(std::size_t ii = 0; ii < ksize; ii++) {
	  //   std::cout << inner_kernel[ii] << ", ";	    
	  // }
	  // std::cout << std::endl;
	  
	  interpolate_1d<KernelT, vec_dims>(inner_kernel, column_start + vec_dims * (inner_int - 1), inner_interp);

	  // std::cout << "inner_interp = ";
	  // for(std::size_t ii = 0; ii < vec_dims; ii++) {
	  //   std::cout << inner_interp[ii] << ", ";	    
	  // }
	  // std::cout << std::endl;

	  scalar_t outer_weight = outer_weights[i1 - (outer_int[0] - 1)][i2 - (outer_int[1] - 1)];
	  // std::cout << "outer_weight = " << outer_weight << std::endl;

	  std::span<scalar_t, vec_dims> output = retval[inner_elem];
	  for(std::size_t cv = 0; cv < vec_dims; cv++) {
	    output[cv] += inner_interp[cv] * outer_weight;
	  }
	  
	  inner_elem++;	  
	}	
      }
    }       
  }
}
