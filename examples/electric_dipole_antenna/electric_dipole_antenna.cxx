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
  CoordVector start_coords = CU::MakeCoordVectorTRZ(0.0, 0.0, -10.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(300.0, 300.0, 30.0);

  // Filter parameters
  scalar_t tp = 5.0;
  unsigned int N = 4;

  // Sampling parameters
  scalar_t os_factor = 30;
  scalar_t r_min = 0.1;

  WFU::CreateElectricDipoleWeightingField(wf_path, start_coords, end_coords, tp, N, r_min, os_factor);
}
