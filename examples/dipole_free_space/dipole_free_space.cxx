#include "Common.hh"
#include "AnalyticGreensFunctionCalculator.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string gf_path = argv[1];

  // Domain of weighting field
  RZCoordVector start_coords{0.0f, -200.0f};
  RZCoordVector end_coords{300.0f, 100.0f};
  scalar_t t_end = 400;

  // Index of refraction
  scalar_t ior = 1.42;
  
  // Filter parameters
  scalar_t filter_t_peak = 3.0;   // lowpass filter peaking time
  unsigned int filter_order = 4;  // lowpass filter order

  // Sampling parameters
  scalar_t os_factor = 15;
  scalar_t r_min = 0.1;

  GreensFunctionCalculator::Analytic::ElectricDipole(gf_path, start_coords, end_coords, t_end, ior, filter_t_peak, filter_order, r_min, os_factor);
}
