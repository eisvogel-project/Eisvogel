#include "Eisvogel/Common.hh"
#include "Eisvogel/GreensFunctionUtils.hh"

namespace GFU = GreensFunctionUtils;

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string gf_path = argv[1];

  // Domain of weighting field
  RZTCoordVector start_coords{0.0f, 0.0f, -50.0f};
  RZTCoordVector end_coords{50.0f, 70.0f, 50.0f};

  // Index of refraction
  scalar_t ior = 1.0;
  
  // Filter parameters
  scalar_t filter_t_peak = 2.0;   // lowpass filter peaking time
  unsigned int filter_order = 4;  // lowpass filter order

  // Sampling parameters
  scalar_t os_factor = 30;
  scalar_t r_min = 0.1;

  GFU::CreateElectricDipoleGreensFunction(gf_path, start_coords, end_coords, ior, filter_t_peak, filter_order, r_min, os_factor);
}
