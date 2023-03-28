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

#include "Serialization.hh"

#include <fstream>

namespace WFU = WeightingFieldUtils;

int main(void) {

  std::string path = "/home/windischhofer/Eisvogel/Eisvogel/signal/build/test.bin";
  std::fstream ofs;
  
  // ofs.open(path, std::ios::out | std::ios::binary);
  ofs.open(path, std::ios::in | std::ios::binary);
  stor::Serializer ser(ofs);

  std::cout << "Building weighting field ..." << std::endl;
  CoordVector start_coords = CoordUtils::MakeCoordVectorTRZ(-10.0, -10.0, -30.0);
  CoordVector end_coords = CoordUtils::MakeCoordVectorTRZ(320.0, 300.0, 30.0);
  scalar_t tp = 5.0;
  unsigned int N = 4;
  // WeightingField wf = WFU::CreateElectricDipoleWeightingField(start_coords, end_coords, tp, N, 70);
  // ser.serialize(wf);
  // ofs.close();

  WeightingField wf = ser.deserialize<WeightingField>();

  SplineInterpolationKernelOrder1 interpolation_kernel;
  // SplineInterpolationKernelOrder3 interpolation_kernel;
  // SincInterpolationKernel interpolation_kernel;

  Integrator integrator(wf, interpolation_kernel);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 1;
  scalar_t tstart = -250, tend = 250;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;
  Trajectory traj({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	}
    );

  std::cout << "Computing signal ..." << std::endl;

  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = 2; cur_t < 3; cur_t += 1) {
    scalar_t cur_signal = integrator.integrate(cur_t, traj);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
