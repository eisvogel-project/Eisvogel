#include <iostream>
#include <cmath>
#include <vector>

#include "Common.hh"
#include "NDArray.hh"
#include "WeightingFieldUtils.hh"
#include "Interpolator.hh"
#include "Integrator.hh"

#include "Trajectory.hh"

#include "Kernels.hh"
#include "IteratorUtils.hh"

#include "MathUtils.hh"

#include "SignalExport.hh"

namespace WFU = WeightingFieldUtils;

int main(void) {

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf = WFU::CreateElectricDipoleWeightingField();
  SplineInterpolationKernelOrder1 interpolation_kernel;
  Integrator integrator(wf, interpolation_kernel);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 1;
  scalar_t tstart = -300, tend = 300;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;
  Trajectory traj({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	}
    );

  std::cout << "Computing signal ..." << std::endl;

  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = -15; cur_t < 15; cur_t += 1) {
    scalar_t cur_signal = integrator.integrate(cur_t, traj);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
