#include <iostream>
#include <cmath>
#include "Common.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"
#include "Antenna.hh"
#include "Geometry.hh"

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  if(argc != 3) {
    std::cout << "Usage: dipole_ice OUTPUT_PATH SCRATCH_PATH" << std::endl;
    return 0;
  }

  std::filesystem::path gf_path = argv[1];
  std::filesystem::path scratch_path = argv[2];
  
  auto eps = []([[maybe_unused]] scalar_t r, scalar_t z) {
    
    scalar_t z_m = z / 3.0;    
    if(z_m > 0.0) {
      return 1.0;
    }

    scalar_t density = 0.0;

    if(z_m > -14.9) {
      density = 0.917 - 0.594 * std::exp(z_m / 30.8);
    }
    else {
      density = 0.917 - 0.367 * std::exp((z_m + 14.9) / 40.5);
    }

    double eps = std::pow(1 + 0.845 * density, 2.0);
    return eps;
  };

  auto impulse_response = [](scalar_t t) {
    unsigned int N = 4; // order of filter
    double tp = 2.0; // peaking time of filter
    if(t <= 0) {
      return 0.0;
    }
    return std::pow(t / tp * N, N) * std::exp(-t / tp * N) / (tp * std::exp(std::lgamma(N)));
  };
  
  CylinderGeometry geom(300, -300, 200, eps);
  InfEDipoleAntenna dipole(0.0, 50.0, -100.0, impulse_response);
  scalar_t t_end = 300;
  
  GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator gfc(geom, dipole, t_end);
  gfc.Calculate(gf_path, scratch_path, scratch_path);
    
  std::cout << "done" << std::endl;
  
  return 0;
}
