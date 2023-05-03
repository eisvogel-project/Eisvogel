#include "Eisvogel/Kernels.hh"
#include <cmath>

std::size_t KeysCubicInterpolationKernel::Support() const {
  return 2;
}

scalar_t KeysCubicInterpolationKernel::operator()(scalar_t arg) const {
  scalar_t abs_arg = std::fabs(arg);

  if(abs_arg < 1.0) {
    return 1.0 + abs_arg * abs_arg * (-2.5 + 1.5 * abs_arg);
  }
  else if(abs_arg < 2.0) {
    return 2.0 + abs_arg * (-4.0 + (2.5 - 0.5 * abs_arg) * abs_arg);
  }
  else {
    return 0.0;
  }
}

scalar_t KeysCubicInterpolationKernel::CDF(int arg) const {

  if(arg <= -2) {
    return 0.0;
  }
  else if(arg == -1) {
    return -0.0416667;
  }
  else if(arg == 0) {
    return 0.5;
  }
  else if(arg == 1) {
    return 1.04167;
  }
  else if(arg >= 2) {
    return 1.0;
  }

  return 1.0;
}
