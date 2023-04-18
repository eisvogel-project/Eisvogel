#include <iostream>
#include <fstream>

#include "Eisvogel/Common.hh"
#include "Eisvogel/Serialization.hh"
#include "Eisvogel/WeightingFieldUtils.hh"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];

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
}
