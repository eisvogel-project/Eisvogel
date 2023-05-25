#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingFieldUtils.hh"

namespace WFU = WeightingFieldUtils;
namespace CU = CoordUtils;

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];

  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-50.0, -10.0, -10.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(50.0, 10.0, 10.0);

  // Filter parameters
  scalar_t tp = 1.0;
  unsigned int N = 6;

  // Sampling parameters
  scalar_t os_factor = 30;
  scalar_t r_min = 0.1;
  
  scalar_t index_of_refraction = 1.3;

  std::cout << "Building weighting field ..." << std::endl;
  WeightingField wf_out = WFU::CreateElectricDipoleWeightingField(start_coords, end_coords, tp, N, r_min, os_factor, index_of_refraction);

  std::cout << "Saving weighting field ..." << std::endl;
  oser.serialize(wf_out);
  ofs.close();
}
