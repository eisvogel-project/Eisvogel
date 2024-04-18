#include "Eisvogel/Common.hh"
#include "Eisvogel/GreensFunctionUtils.hh"

/*
Script to generate a Green's function that covers the time and space
required to perform the shower simulation.
*/

namespace GFU = GreensFunctionUtils;

int main(int argc, char* argv[]) {

  std::string gf_path;
  if (argc < 2) {
    gf_path = "greens_function";
  } else {
    gf_path = argv[1];
  }

  // Domain of Green's function
  RZCoordVector start_coords{0.0f, -551.0f};
  RZCoordVector end_coords{700.0f, -549.0f};
  scalar_t t_end = 1350.0;
  
  // Filter parameters
  scalar_t filter_t_peak = 1.0;
  unsigned int filter_order = 6;

  // Sampling parameters
  scalar_t os_factor = 5;
  scalar_t r_min = 0.1;
  
  scalar_t index_of_refraction = 1.78;

  std::cout << "Building Green's function ..." << std::endl;

  GFU::CreateElectricDipoleGreensFunction(gf_path, start_coords, end_coords, t_end, index_of_refraction, filter_t_peak, filter_order, r_min, os_factor);
}
