#include <iostream>
#include <fstream>
#include <algorithm>

#include "Eisvogel/Common.hh"
#include "Eisvogel/Integrator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/Kernels.hh"
#include "Eisvogel/SignalExport.hh"
#include "Eisvogel/Serialization.hh"
#include "Eisvogel/WeightingFieldUtils.hh"
#include "shower_creator.h"
#include "shower_1D.h"
#include "constants.h"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(void) {
    
  std::string wf_path = "electric_dipole_wf.bin";

  std::fstream ofs;
  ofs.open(wf_path, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);

  //std::array<float, 3> shower_vertex = {-76, .1, -63};
  std::array<float, 3> shower_vertex = {-76, .1, -63};
  std::cout << "Building Shower \n";

  showers::ShowerCreator shower_creator("/home/welling/RadioNeutrino/scripts/Eisvogel/extern/shower_profile/shower_file");
  showers::Shower1D shower = shower_creator.create_shower(
          shower_vertex,
          2.0e+19,
          1.5708,
          0,
          1
          );

  std::vector<double> t;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> ce;
  shower.get_shower(1., &t, &x, &y, &z, &ce);
  double r_max = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double r_min = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double z_max = z[0];
  double z_min = z[0];
  double r;
  for (int i=0; i<t.size(); i++) {
      r = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
      r_min = std::min(r_min, r);
      r_max = std::max(r_max, r);
      z_min = std::min(z_min, z[i]);
      z_max = std::max(z_max, z[i]);
  }
  r_max /= constants::c;
  r_min /= constants::c;
  z_max /= constants::c;
  z_min /= constants::c;
  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-500, -100, z_min - 200);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(1000.0, r_max + 800, z_max + 200);

  // Filter parameters
  scalar_t tp = 5.0;
  unsigned int N = 4;

  // Sampling parameters
  scalar_t os_factor = 7;
  scalar_t R_min = 0.1;
  scalar_t refractive_index = 1.3;

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf_out = WFU::CreateElectricDipoleWeightingField(
          start_coords,
          end_coords,
          tp,
          N,
          R_min,
          os_factor,
          refractive_index
          );

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
  
  
  std::vector<scalar_t> signal_times, signal_values;
  std::cout << "Integrating \n";
  for(scalar_t cur_t = 0; cur_t < 800; cur_t += 1) {
    scalar_t cur_signal = integrator.integrate(cur_t, shower, 1);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }
  ExportSignal(signal_times, signal_values, "./2e19HAD.csv");
  return 0;
}
