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

    // weight for sample at `i_int - 1`
    // *(out+0) = i_frac * (-0.5 + (1.0 - 0.5 * i_frac) * i_frac);
    scalar_t p1_1 = -0.5 * i_frac;
    scalar_t p1_2 = i_frac - 1.0;
    scalar_t p1_3 = p1_1 * p1_2;
    *(out+0) = p1_3 * p1_2;

    // weight for sample at `i_int`
    // *(out+1) = 1.0 + i_frac * i_frac * (-2.5 + 1.5 * i_frac);
    constexpr scalar_t c2_1 = (1.0 - std::sqrt(7.0)) / 3.0;
    scalar_t p2_1 = i_frac - c2_1;
    constexpr scalar_t c2_2 = (1.0 + std::sqrt(7.0)) / 3.0;
    scalar_t p2_2 = i_frac - c2_2;
    scalar_t p2_3 = i_frac - 1.0;
    scalar_t p2_4 = p2_2 * p2_3;
    constexpr scalar_t c2_5 = 3.0 / 2.0;
    scalar_t p2_5 = c2_5 * p2_1;
    *(out+1) = p2_5 * p2_4;

    // weight for sample at `i_int + 1`
    // *(out+2) = i_frac * (0.5 + (2.0 - 1.5 * i_frac) * i_frac);
    constexpr scalar_t c3_1 = 1.0 / 3.0 * (2.0 - std::sqrt(7.0));
    scalar_t p3_1 = i_frac - c3_1;
    constexpr scalar_t c3_2 = 1.0 / 3.0 * (2.0 + std::sqrt(7.0));
    scalar_t p3_2 = i_frac - c3_2;
    scalar_t p3_3 = p3_1 * p3_2;
    constexpr scalar_t c3_4 = -3.0 / 2.0;
    scalar_t p3_4 = i_frac * c3_4;        
    *(out+2) = p3_3 * p3_4;

    // weight for sample at `i_int + 2`
    // *(out+3) = (-0.5 + 0.5 * i_frac) * i_frac * i_frac;
    scalar_t p4_1 = i_frac - 1.0;
    scalar_t p4_2 = i_frac * i_frac;
    constexpr scalar_t c4_3 = 1.0 / 2.0;
    scalar_t p4_3 = p4_1 * c4_3;
    *(out+3) = p4_3 * p4_2;
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
		   scalar_t inner_ind_start, scalar_t inner_ind_delta, std::size_t num_samples) {

    // Ensure that the output array can hold all interpolated points
    retval.resize(num_samples);
    
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

    constexpr std::size_t ind_offset = KernelT::support - 1;
    
    for(std::size_t i1 = outer_int[0] - ind_offset; i1 < outer_int[0] - ind_offset + ksize; i1++) {
      for(std::size_t i2 = outer_int[1] - ind_offset; i2 < outer_int[1] - ind_offset + ksize; i2++) {

	Vector<std::size_t, dims> column_ind(0);
	column_ind[0] = i1;
	column_ind[1] = i2;
	scalar_t* column_start = &(*(arr.m_data -> begin())) + arr.ComputeFlatInd(column_ind);

	scalar_t outer_weight = outer_weights[i1 - (outer_int[0] - ind_offset)][i2 - (outer_int[1] - ind_offset)];
	// std::cout << "outer_weight = " << outer_weight << std::endl;
	
	for(std::size_t cur_sample_ind = 0; cur_sample_ind < num_samples; cur_sample_ind++) {
	  //for(scalar_t cur_inner_ind = inner_ind_start; cur_inner_ind < inner_ind_end; cur_inner_ind += inner_ind_delta) {

	  scalar_t cur_inner_ind = inner_ind_start + cur_sample_ind * inner_ind_delta;
	  // std::cout << "on inner_elem = " << inner_elem << std::endl;
	  
	  scalar_t tmp = 0.0;
	  scalar_t inner_frac = std::modf(cur_inner_ind, &tmp);
	  std::size_t inner_int = (std::size_t)tmp;

	  // std::cout << "inner_frac = " << inner_frac << std::endl;

	  // at the moment, this recomputes the kernel for every inner interpolation -> need to batch this
	  scalar_t inner_interp[vec_dims] = {0};
	  scalar_t inner_kernel[ksize] = {0};
	  KernelT::kern(inner_frac, inner_kernel);

	  // std::cout << "inner_kernel = ";
	  // for(std::size_t ii = 0; ii < ksize; ii++) {
	  //   std::cout << inner_kernel[ii] << ", ";	    
	  // }
	  // std::cout << std::endl;
	  
	  interpolate_1d<KernelT, vec_dims>(inner_kernel, column_start + vec_dims * (inner_int - ind_offset), inner_interp);

	  // std::cout << "inner_interp = ";
	  // for(std::size_t ii = 0; ii < vec_dims; ii++) {
	  //   std::cout << inner_interp[ii] << ", ";	    
	  // }
	  // std::cout << std::endl;

	  std::span<scalar_t, vec_dims> output = retval[cur_sample_ind];
	  for(std::size_t cv = 0; cv < vec_dims; cv++) {
	    output[cv] += inner_interp[cv] * outer_weight;
	  }	 
	}	
      }
    }       
  }
}
