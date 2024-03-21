#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingFieldUtilsOld.hh"

namespace WFU = WeightingFieldUtilsOld;
namespace CU = CoordUtils;
/*
Script to generate a weighting field that covers the time and space
required to perform the shower simulation.
*/


int main(int argc, char* argv[]) {

  std::string wf_path;
  if (argc < 2) {
    wf_path = "weighting_field";
  } else {
    wf_path = argv[1];
  }

  // Domain of weighting field
  CoordVector start_coords = CU::MakeCoordVectorTRZ(-500.0, -10.0, -551.0);
  CoordVector end_coords = CU::MakeCoordVectorTRZ(1350.0, 700.0, -549.0);

  // Filter parameters
  scalar_t tp = 1.0;
  unsigned int N = 6;

  // Sampling parameters
  scalar_t os_factor = 5;
  scalar_t r_min = 0.1;
  
  scalar_t index_of_refraction = 1.78;

  std::cout << "Building weighting field ..." << std::endl;

  WFU::CreateElectricDipoleWeightingField(wf_path, start_coords, end_coords, tp, N, r_min, os_factor, index_of_refraction);
}
