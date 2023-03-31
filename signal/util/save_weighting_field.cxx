#include <iostream>
#include <fstream>
#include "Common.hh"
#include "Serialization.hh"
#include "CoordUtils.hh"
#include "WeightingFieldUtils.hh"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(void) {

  std::string path = "/home/windischhofer/Eisvogel/Eisvogel/signal/build/electric_dipole_wf.bin";

  std::fstream ofs;
  ofs.open(path, std::ios::out | std::ios::binary);  
  stor::Serializer ser(ofs);

  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-20.0, -10.0, -40.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(320.0, 300.0, 40.0);

  // Filter parameters
  scalar_t tp = 5.0;
  unsigned int N = 4;

  // Sampling parameters
  scalar_t os_factor = 4;
  scalar_t r_min = 0.1;

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf = WFU::CreateElectricDipoleWeightingField(start_coords, end_coords, tp, N, r_min, os_factor);

  std::cout << "Saving weighting field ..." << std::endl;
  ser.serialize(wf);
  ofs.close();

  return 0;
}
