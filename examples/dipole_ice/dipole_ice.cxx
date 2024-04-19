#include <iostream>
#include <cmath>
#include "Eisvogel/Common.hh"
#include "Eisvogel/CylindricalGreensFunctionCalculator.hh"
#include "Eisvogel/Antenna.hh"
#include "Eisvogel/Geometry.hh"

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  if(argc < 2) {
    throw std::runtime_error("Error: need to pass path to output directory!");
  }

  std::string gf_path = argv[1];
  
  auto eps = [](scalar_t r, scalar_t z) {
    
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

  // CylinderGeometry geom(20.0, -15.0, 15.0, eps);
  // InfEDipoleAntenna dipole(0.0, 10.0, 0.0, impulse_response);
  // scalar_t t_end = 25.0;

  // CylinderGeometry geom(100.0, -100.0, 100.0, eps);
  // InfEDipoleAntenna dipole(0.0, 10.0, 0.0, impulse_response);
  // scalar_t t_end = 25.0;
  
  CylinderGeometry geom(400, -400, 200, eps);
  InfEDipoleAntenna dipole(0.0, 10.0, -100.0, impulse_response);
  scalar_t t_end = 500;
  
  CylindricalGreensFunctionCalculator gfc(geom, dipole, t_end);
  gfc.Calculate(gf_path, "/scratch/midway3/windischhofer/", "/scratch/midway3/windischhofer/");
    
  std::cout << "done" << std::endl;
  
  return 0;
}
