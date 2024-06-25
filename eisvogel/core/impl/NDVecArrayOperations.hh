#pragma once

#include "Vector.hh"
#include "NDVecArray.hh"

namespace Downsampling {

  template <std::size_t dims>
  Vector<std::size_t, dims> get_downsampled_shape(
      const Vector<std::size_t, dims>& start_ind,
      const Vector<std::size_t, dims>& shape,
      const Vector<std::size_t, dims>& downsampling);

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void downsample(const NDVecArray<T, dims, vec_dims>& to_downsample,
                  const Vector<std::size_t, dims> start_ind,
                  const Vector<std::size_t, dims> downsampling,
                  NDVecArray<T, dims, vec_dims>& downsampled,
                  Vector<std::size_t, dims>& downsampled_start_ind);
} // namespace Downsampling

#include "NDVecArrayOperations.hxx"
