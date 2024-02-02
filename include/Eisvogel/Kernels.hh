#ifndef __KERNELS_HH
#define __KERNELS_HH

#include <limits>
#include <cmath>
#include "Common.hh"

struct Kernel {
  
  virtual std::size_t Support() const = 0;
  virtual scalar_t operator()(scalar_t arg) const = 0;
  virtual scalar_t CDF(int arg) const = 0;
};

struct KeysCubicInterpolationKernelNew {

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
};

struct KeysCubicInterpolationKernel : public Kernel {

  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(int arg) const;
};

#endif
