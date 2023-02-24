#ifndef __MATH_UTILS_HH
#define __MATH_UTILS_HH

#include "Common.hh"

namespace MathUtils {

  scalar_t incomplete_gamma(scalar_t a, scalar_t z, scalar_t rel_tol = 1e-6);

  scalar_t incomplete_gamma_series_expansion(scalar_t a, scalar_t z, scalar_t rel_tol);
  scalar_t incomplete_gamma_continued_fraction(scalar_t a, scalar_t z, scalar_t rel_tol);

};

#endif
