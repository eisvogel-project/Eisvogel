#include <iostream>
#include <fstream>

#include "Common.hh"
#include "Integrator.hh"
#include "Current0D.hh"
#include "Kernels.hh"
#include "SignalExport.hh"
#include "Serialization.hh"
#include "WeightingFieldUtils.hh"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(void) {

  std::string wf_path = "electric_dipole_wf.bin";

  std::fstream ofs;
  ofs.open(wf_path, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);

  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-20.0, -10.0, -40.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(320.0, 500.0, 40.0);

  // Filter parameters
  scalar_t tp = 5.0;
  unsigned int N = 4;

  // Sampling parameters
  scalar_t os_factor = 7;
  scalar_t r_min = 0.1;

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf_out = WFU::CreateElectricDipoleWeightingField(start_coords, end_coords, tp, N, r_min, os_factor);

  std::cout << "Saving weighting field ..." << std::endl;
  oser.serialize(wf_out);
  ofs.close();

  std::fstream ifs; 
  ifs.open(wf_path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);

  std::cout << "Loading weighting field ..." << std::endl;
  WeightingField wf_in = iser.deserialize<WeightingField>();

  KeysCubicInterpolationKernel interpolation_kernel;
  Integrator integrator(wf_in, interpolation_kernel);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = -250, tend = 250;
  scalar_t charge = 1;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;
  Current0D curr({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	},
    {charge}
    );

  std::cout << "Computing signal ..." << std::endl;

  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = -10; cur_t < 20; cur_t += 1) {
    scalar_t cur_signal = integrator.integrate(cur_t, curr);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
