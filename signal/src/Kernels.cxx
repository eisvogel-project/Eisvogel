#include "Kernels.hh"
#include <cmath>

std::size_t SplineInterpolationKernelOrder1::Support() const {
  return 1;
}

// Can probably make this more efficient with some sort of switch/case on the integer part of the argument
scalar_t SplineInterpolationKernelOrder1::operator()(scalar_t arg) const {
  scalar_t abs_arg = std::fabs(arg);

  if(InRange(0, 1, abs_arg)) {
    return 1.0 - abs_arg;
  }
  else if(InRange(1, Inf, abs_arg)) {
    return 0.0;
  }
  
  throw;
}

std::size_t SplineInterpolationKernelOrder3::Support() const {
  return 1;
}

scalar_t SplineInterpolationKernelOrder3::operator()(scalar_t arg) const {
  scalar_t abs_arg = std::fabs(arg);

  return 0.0;

  throw;
}
