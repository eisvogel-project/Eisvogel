#include <iostream>
#include "Common.hh"
#include "NDArray.hh"
#include "WeightingFieldUtils.hh"
#include "Interpolator.hh"

#include "Kernels.hh"

#include <cmath>
#include <tuple>

namespace WFU = WeightingFieldUtils;

template <typename... FracInds>
scalar_t Interpolate(FracInds... frac_inds) {

  std::size_t support = 4;

  std::array<scalar_t, sizeof... (FracInds)> central_inds{{frac_inds...}};

  for(auto cur: central_inds) {
    std::cout << cur << " ";
  }
  std::cout << std::endl;

  return 0.0;
}

int main(void) {

  unsigned int len_t = 10;

  int n = 10;

  DenseWeightingField wf = WFU::CreateElectricDipoleWeightingField();

  std::cout << wf.E_r(0,0,0) << std::endl;
  std::cout << wf.E_z(0,0,0) << std::endl;
  std::cout << wf.E_phi(0,0,0) << std::endl;

  DenseNDArray<scalar_t, 3> testarr({10, 10, 10});

  Interpolator<DenseNDArray, scalar_t, 3> itpl(testarr);

  std::cout << itpl.Interpolate(2, 2, 2) << std::endl;

  Kernels::SplineInterpolationKernelOrder1 kernel;

  std::cout << kernel.Support() << std::endl;
  std::cout << kernel(0.78) << std::endl;
  std::cout << kernel(1.78) << std::endl;

  std::cout << "-----" << std::endl;

  std::cout << Interpolate(1.234f, 2.53f) << std::endl;

  return 0;
}
