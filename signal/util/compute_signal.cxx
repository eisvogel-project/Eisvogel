#include <iostream>
#include <fstream>

#include "Common.hh"
#include "Integrator.hh"
#include "Current0D.hh"
#include "Kernels.hh"
#include "SignalExport.hh"
#include "Serialization.hh"

int main(void) {

  std::string path = "/home/windischhofer/Eisvogel/Eisvogel/signal/build/electric_dipole_wf.bin";
  std::fstream ofs;
  
  ofs.open(path, std::ios::in | std::ios::binary);
  stor::Serializer ser(ofs);

  std::cout << "Loading weighting field ..." << std::endl;
  WeightingField wf = ser.deserialize<WeightingField>();

  // SplineInterpolationKernelOrder1 interpolation_kernel;
  // SplineInterpolationKernelOrder3 interpolation_kernel;
  KeysCubicInterpolationKernel interpolation_kernel;
  Integrator integrator(wf, interpolation_kernel);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 20;
  scalar_t tstart = -250, tend = 250;
  scalar_t charge = 1;
  scalar_t beta = 1.2;

  std::cout << "Building trajectory ..." << std::endl;
  Current0D curr({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	},
    {charge}
    );

  std::cout << "Computing signal ..." << std::endl;

  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = 0; cur_t < 40; cur_t += 1) {
    scalar_t cur_signal = integrator.integrate(cur_t, curr);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
