#include <iostream>
#include <cmath>
#include "Eisvogel/Common.hh"
#include "Eisvogel/CylindricalWeightingFieldCalculator.hh"
#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

unsigned int fact(unsigned arg) {
  unsigned int retval = 1;

  for(unsigned int cur = 1; cur <= arg; cur++) {
    retval *= cur;
  }
  
  return retval;
}

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  if(argc < 2) {
    throw std::runtime_error("Error: need to pass path to output directory!");
  }

  std::string wf_path = argv[1];
  
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

  // CylinderGeometry geom(20, -20, 20, eps);
  // InfEDipoleAntenna dipole(0.0, 10.0, -5.0, impulse_response);
  
  CylinderGeometry geom(600, -300, 100, eps);
  InfEDipoleAntenna dipole(0.0, 10.0, -100.0, impulse_response);

  scalar_t t_end = 25;
  CylindricalWeightingFieldCalculator wfc(geom, dipole, t_end);
  wfc.Calculate(wf_path);
  
  return 0;
}
