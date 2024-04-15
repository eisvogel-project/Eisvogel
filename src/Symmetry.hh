#pragma once

#include "DistributedNDVecArray.hh"

namespace SpatialSymmetry {

  template <typename T, std::size_t vec_dims>
  class Cylindrical {

  private:

    static constexpr std::size_t dims = 3;        // A cylindrically-symmetric Green's function is indexed as (R, Z, T) for m = 0
    
  public:

    using darr_t = DistributedNDVecArray<NDVecArray, T, dims, vec_dims>;
    using chunk_t = NDVecArray<T, dims, vec_dims>;
    using view_t = typename darr_t::view_t;
    
    static void boundary_evaluator(darr_t& darr, const RZTSignedIndexVector& ind, view_t elem);
    
  };
}

#include "Symmetry.hxx"
