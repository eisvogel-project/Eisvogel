#pragma once

#include "NDVecArray.hh"
#include "Vector.hh"
#include "Common.hh"

namespace Interpolation::Kernel {

  struct Keys {
    static constexpr std::size_t support = 2;
    static void kern(scalar_t i_frac, scalar_t* out);
  };

  struct Linear {
    static constexpr std::size_t support = 1;
    static void kern(scalar_t i_frac, scalar_t* out);
  };
} // namespace Interpolation::Kernel

namespace Interpolation {

  template <class KernelT, typename T, std::size_t dims, std::size_t vec_dims>
  void interpolate(const NDVecArray<T, dims, vec_dims>& arr,
                   NDVecArray<T, 1, vec_dims>& retval,
                   const Vector<scalar_t, dims - 1>& outer_inds,
                   scalar_t inner_ind_start, scalar_t inner_ind_delta,
                   std::size_t num_samples);

}

#include "Interpolation.hxx"
