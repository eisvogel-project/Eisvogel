#ifndef __KERNELS_HH
#define __KERNELS_HH

#include <limits>
#include <cmath>
#include "Common.hh"

struct KeysCubicInterpolationKernel {

  const static int Support = 2;
  
  static scalar_t Weight(scalar_t arg) {
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

  static scalar_t CDF(int arg) {

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
};

using DefaultKernel = KeysCubicInterpolationKernel;

#endif
