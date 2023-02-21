#include <iostream>
#include <cmath>
#include <tuple>

#include "Common.hh"
#include "NDArray.hh"
#include "WeightingFieldUtils.hh"
#include "Interpolator.hh"
#include "Integrator.hh"

#include "Trajectory.hh"

#include "Kernels.hh"
#include "IteratorUtils.hh"

namespace WFU = WeightingFieldUtils;

int main(void) {

  unsigned int len_t = 10;

  int n = 10;

  DenseNDArray<scalar_t, 1> vector1({1, 2, 3});
  DenseNDArray<scalar_t, 1> vector2({1, 2, 3});

  DenseNDArray<scalar_t, 1> vector3 = (vector1 + 2 + vector2) / 3;

  for(auto cur: vector3) {
    std::cout << cur << " ";
  }
  std::cout << std::endl;

  DenseNDArray<scalar_t, 3> testarr({10, 10, 10}, 2.0);
  testarr(2, 3, 1) = 18;
  IndexVector inds({2, 3, 1});

  std::cout << "HHH " << testarr(inds) << " HHH" << std::endl;

  SplineInterpolationKernelOrder1 kernel;

  std::cout << kernel.Support() << std::endl;
  std::cout << kernel(0.78) << std::endl;
  std::cout << kernel(1.78) << std::endl;

  std::cout << "-----" << std::endl;

  DenseNDArray<scalar_t, 1> testarr1d({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  Interpolator<DenseNDArray, scalar_t, 1> itpl(testarr1d, kernel);

  std::cout << itpl.Interpolate(0.5) << std::endl;

  std::cout << "-----" << std::endl;

  WeightingField wf = WFU::CreateElectricDipoleWeightingField();
  SplineInterpolationKernelOrder1 interpolation_kernel;

  Integrator integrator(wf, interpolation_kernel);

  // ======

  scalar_t b = 10;
  scalar_t tstart = -300, tend = 300;
  scalar_t beta = 0.9;

  Trajectory traj({
      CoordUtils::MakeCoordVectorTXYZ(tstart, 2, beta * tstart, b),
      CoordUtils::MakeCoordVectorTXYZ(tend, 2, beta * tend, b)
	}
    );

  integrator.integrate(0.0, traj);

  // =======

  return 0;
}
