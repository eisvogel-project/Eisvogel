#include <iostream>
#include <fstream>

#include "Eisvogel/Common.hh"
#include "Eisvogel/Integrator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/Kernels.hh"
#include "Eisvogel/SignalExport.hh"
#include "Eisvogel/Serialization.hh"
#include "Eisvogel/WeightingFieldUtils.hh"
#include "shower_creator.h"
#include "shower_1D.h"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(void) {
    
  std::string wf_path = "electric_dipole_wf.bin";

  std::fstream ofs;
  ofs.open(wf_path, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);

  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-1000.0, -500.0, -50.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(1000.0, 1000.0, 50.0);

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
  scalar_t tstart = 0, tend = 20;
  scalar_t charge = 1;
  scalar_t beta = 0.9;
  
  std::array<float, 3> shower_vertex = {0, 0, -10};
  std::cout << "Building Shower \n";

  showers::ShowerCreator shower_creator("/home/welling/RadioNeutrino/scripts/Eisvogel/extern/shower_profile/shower_file");
  showers::Shower1D shower = shower_creator.create_shower(
          shower_vertex,
          5.0e+17,
          1.5708,
          0,
          true
          );

  std::vector<scalar_t> signal_times, signal_values;
  std::cout << "Interating \n";
  for(scalar_t cur_t = 0; cur_t < 10; cur_t += 1) {
      std::cout << "Calculating t=" << cur_t << "\n";
    scalar_t cur_signal = integrator.integrate(cur_t, shower, 1);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }
  ExportSignal(signal_times, signal_values, "./test_signal.csv");
  return 0;
}