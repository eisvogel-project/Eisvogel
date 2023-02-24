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

  // Algorithm from https://doi.org/10.2307/2347328
  scalar_t incomplete_gamma_continued_fraction(scalar_t a, scalar_t z, scalar_t rel_tol) {

    scalar_t An2 = 0;
    scalar_t An1 = 1;

    scalar_t Bn2 = 1;
    scalar_t Bn1 = z - a + 1;

    scalar_t an = a - 1;
    scalar_t bn = Bn1 + 2;

    std::size_t n = 1;
    scalar_t contfrac = 0.0;

    while(true) {
      
      scalar_t An = bn * An1 + an * An2;
      scalar_t Bn = bn * Bn1 + an * Bn2;

      An2 = An1;
      An1 = An;

      Bn2 = Bn1;
      Bn1 = Bn;

      scalar_t next_contfrac = An / Bn;
      scalar_t rel_err = std::fabs(next_contfrac - contfrac) / next_contfrac;
      contfrac = next_contfrac;

      if(rel_err < rel_tol) {
	break;
      }      

      an += (a - 2 * n - 1);
      bn += 2;
      n += 1;
    }

    return std::exp(a * std::log(z) - z) * contfrac;
  }

  scalar_t incomplete_gamma(scalar_t a, scalar_t z, scalar_t rel_tol) {    

    if((a < 0) || (z < 0)) {
      throw;
    }

    if((z <= 1) || (z < a)) {
      return incomplete_gamma_series_expansion(a, z, rel_tol);
    }
    else {
      return incomplete_gamma_continued_fraction(a, z, rel_tol);
    }    
  }
}
