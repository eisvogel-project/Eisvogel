#pragma once

#include <cstddef>
#include "Common.hh"

namespace MathUtils {

  bool is_integer(scalar_t k) {
    return std::floor(k) == k;
  }
  
  constexpr std::size_t IntDivCeil(std::size_t x, std::size_t y) {
    // calculates ceil(x/y)
    return (x + y - 1) / y;
  }
  
  scalar_t incomplete_gamma(scalar_t a, scalar_t z, scalar_t rel_tol = 1e-6);

  scalar_t incomplete_gamma_series_expansion(scalar_t a, scalar_t z, scalar_t rel_tol);
  scalar_t incomplete_gamma_continued_fraction(scalar_t a, scalar_t z, scalar_t rel_tol);

};

#include "MathUtils.hxx"
