#include "Kernels.hh"
#include <cmath>

std::size_t SincInterpolationKernel::Support() const {
  return 2; // some cutoff
}

scalar_t SincInterpolationKernel::operator()(scalar_t arg) const {
  scalar_t scaled_arg = M_PI * std::fabs(arg);
  return std::sin(scaled_arg) / scaled_arg;
}

scalar_t SincInterpolationKernel::CDF(scalar_t arg) const {
  return 1.0;
}

std::size_t SplineInterpolationKernelOrder1::Support() const {
  return 1;
}

// Can probably make this more efficient with some sort of switch/case on the integer part of the argument
scalar_t SplineInterpolationKernelOrder1::operator()(scalar_t arg) const {
  scalar_t abs_arg = std::fabs(arg);

  if(InRangeInclusive(0, 1, abs_arg)) {
    return 1.0 - abs_arg;
  }
  else if(InRangeInclusive(1, Inf, abs_arg)) {
    return 0.0;
  }
  
  throw;
}

// TOOD: can simplify this, as will only ever need to evaluate it at integer values (in terms of number of sampling steps)
scalar_t SplineInterpolationKernelOrder1::CDF(scalar_t arg) const {

  if(arg > 0) {
    return 1.0 - CDF(-arg);
  }
  
  if(InRangeInclusive(-1, 0, arg)) {
    return std::pow(arg, 2) / 2.0 + arg + 0.5;
  }

  return 0.0;
}

std::size_t SplineInterpolationKernelOrder3::Support() const {
  return 4;
}

scalar_t SplineInterpolationKernelOrder3::operator()(scalar_t arg) const {
  scalar_t abs_arg = std::fabs(arg);

  if(InRangeInclusive(0, 1, abs_arg)) {
    return 1.0 + std::pow(abs_arg, 2) * (-2.196152422706632 + 1.196152422706632 * abs_arg);
  }
  else if(InRangeInclusive(1, 2, abs_arg)) {
    return 2.7846096908265263 + abs_arg * (-5.353829072479582 + (3.1576766497729514 - 0.5884572681198961 * abs_arg) * abs_arg);
  }
  else if(InRangeInclusive(2, 3, abs_arg)) {
    return -3.1844616523161946 + abs_arg * (3.599777942234521 + (-1.3191268575841093 + 0.1576766497729487 * abs_arg) * abs_arg);
  }
  else if(InRangeInclusive(3, 4, abs_arg)) {
    return 2.2135398277940794 + abs_arg * (-1.798223537875998 + (0.4802069691194646 - 0.04224933097189876 * abs_arg) * abs_arg);
  }
  else if(InRangeInclusive(4, Inf, abs_arg)) {
    return 0.0;
  }

  throw;
}

scalar_t SplineInterpolationKernelOrder3::CDF(scalar_t arg) const {
  return 1.0;
}
