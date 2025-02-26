#pragma once

#include "Vector.hh"

namespace ChunkUtils {

  template <typename T, std::size_t dims>
  bool contains(const Vector<T, dims>& chunk_start_ind, const Vector<T, dims>& chunk_shape, const Vector<T, dims>& ind) {
    for(std::size_t i = 0; i < dims; i++) {
      if((ind[i] - chunk_start_ind[i]) >= chunk_shape[i]) {
	return false;
      }
    }
    return true;
  }
  
  template <typename T, std::size_t dims>
  bool overlaps(const Vector<T, dims>& chunk_a_start_ind, const Vector<T, dims>& chunk_a_shape,
		const Vector<T, dims>& chunk_b_start_ind, const Vector<T, dims>& chunk_b_shape) {

    for(std::size_t i = 0; i < dims; i++) {

      if(chunk_a_start_ind[i] >= chunk_b_start_ind[i] + chunk_b_shape[i]) {
	return false;
      }

      if(chunk_a_start_ind[i] + chunk_a_shape[i] <= chunk_b_start_ind[i]) {
	return false;
      }
    }

    return true;
  }
  
  template <class IndexT, class ShapeT>
  void get_chunk_overlap(const IndexT& chunk_a_start_ind, const ShapeT& chunk_a_shape,
			 const IndexT& chunk_b_start_ind, const ShapeT& chunk_b_shape,
			 IndexT& overlap_start_ind, ShapeT& overlap_shape) {
    overlap_start_ind = VectorUtils::min(chunk_a_start_ind + chunk_a_shape,
					 VectorUtils::max(chunk_a_start_ind, chunk_b_start_ind));
    IndexT overlap_end_ind = VectorUtils::min(chunk_a_start_ind + chunk_a_shape,
					      chunk_b_start_ind + chunk_b_shape);
    overlap_shape = overlap_end_ind - overlap_start_ind;
  }
};
