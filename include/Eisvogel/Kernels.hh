#ifndef __KERNELS_HH
#define __KERNELS_HH

#include <limits>
#include "Common.hh"

class Kernel {
  
public:
  virtual std::size_t Support() const = 0;
  virtual scalar_t operator()(scalar_t arg) const = 0;
  virtual scalar_t CDF(int arg) const = 0;
};

class SincInterpolationKernel : public Kernel {

public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(int arg) const;
};

class SplineInterpolationKernelOrder1 : public Kernel {
  
public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(int arg) const;
};

class SplineInterpolationKernelOrder3 : public Kernel {

public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(int arg) const;
};

class KeysCubicInterpolationKernel : public Kernel {

public:
  std::size_t Support() const;
  scalar_t operator()(scalar_t arg) const;
  scalar_t CDF(int arg) const;
};

#endif
