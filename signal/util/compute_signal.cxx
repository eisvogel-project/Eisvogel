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

#include "MathUtils.hh"

namespace WFU = WeightingFieldUtils;

int main(void) {

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf = WFU::CreateElectricDipoleWeightingField();
  SplineInterpolationKernelOrder1 interpolation_kernel;
  Integrator integrator(wf, interpolation_kernel);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = -300, tend = 300;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;
  Trajectory traj({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	}
    );

  std::cout << "Computing signal ..." << std::endl;
  scalar_t signal = integrator.integrate(0.0, traj);

  return 0;
}
