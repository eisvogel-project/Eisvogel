#include "MathUtils.hh"
#include <cmath>

#include <iostream>

namespace MathUtils {

  // Algorithm from https://doi.org/10.2307/2346431
  scalar_t incomplete_gamma_series_expansion(scalar_t a, scalar_t z, scalar_t rel_tol) {

    scalar_t C = 1.0;
    scalar_t accum = C; 
    scalar_t den = a;

    do {
      den += 1;
      C *= z / den;
      accum += C;
    }
    while(C / accum > rel_tol);

    return std::exp(std::lgamma(a)) * (1 - std::exp(a * std::log(z) - std::lgamma(a + 1) - z) * accum);
  }

  scalar_t incomplete_gamma_continued_fraction(scalar_t a, scalar_t z, scalar_t rel_tol) {
    return 0.0;
  }

  scalar_t incomplete_gamma(scalar_t a, scalar_t z, scalar_t rel_tol) {    

    if((z <= 1) || (z < a)) {
      return incomplete_gamma_series_expansion(a, z, rel_tol);
    }
    else {
      return incomplete_gamma_continued_fraction(a, z, rel_tol);
    }    
  }
}
