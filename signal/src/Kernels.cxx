#include "Kernels.hh"

std::size_t SplineInterpolationKernelOrder1::Support() const {
  return 1;
}

scalar_t SplineInterpolationKernelOrder1::operator()(scalar_t arg) const {
  if(InRange(NegInf, -1, arg)) {
    return 0.0;
  }
  else if(InRange(-1, 0, arg)) {
    return arg + 1.0;
  }
  else if(InRange(0, 1, arg)) {
    return 1.0 - arg;
  }
  else if(InRange(1, Inf, arg)) {
    return 0.0;
  }
  
  throw;
}
