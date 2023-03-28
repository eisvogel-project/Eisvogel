#ifndef __KERNELS_HH
#define __KERNELS_HH

#include <limits>
#include "Common.hh"

inline bool InRangeInclusive(scalar_t min, scalar_t max, scalar_t& val) {
  return (min <= val) && (val <= max);
}

class Kernel {
  
public:
  virtual std::size_t Support() const = 0;
  virtual scalar_t operator()(scalar_t arg) const = 0;
  virtual scalar_t CDF(scalar_t arg) const = 0;
  
protected:
  static constexpr scalar_t Inf = std::numeric_limits<scalar_t>::infinity();
  static constexpr scalar_t NegInf = Inf * (-1);    
};

class SincInterpolationKernel : public Kernel {

public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(scalar_t arg) const;
};

class SplineInterpolationKernelOrder1 : public Kernel {
  
public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(scalar_t arg) const;
};

class SplineInterpolationKernelOrder3 : public Kernel {

public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(scalar_t arg) const;
};

#endif
