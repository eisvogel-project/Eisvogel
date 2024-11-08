#pragma once

#include "DistributedNDVecArray.hh"

namespace SpatialSymmetry {

  // Describes cylindrical symmetry for m = 0
  template <typename T>
  struct Cylindrical {

    static constexpr std::size_t dims = 3;        // A cylindrically-symmetric Green's function is indexed as (R, Z, T) for m = 0 ...
    static constexpr std::size_t vec_dims = 2;    // .. and stores E_r / E_z at each location

    using darr_t = DistributedNDVecArray<NDVecArray, T, dims, vec_dims>;
    using chunk_t = NDVecArray<T, dims, vec_dims>;
    using view_t = typename darr_t::view_t;

    static void boundary_evaluator(chunk_t& chunk, const RZTSignedIndexVector& chunk_start_ind, const RZTSignedIndexVector& chunk_end_ind,
				   const RZTSignedIndexVector& elem_ind, view_t elem);
  };
}

#include "Symmetry.hxx"
