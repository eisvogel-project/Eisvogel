#include <iostream>
#include <cmath>
#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingFieldCalculator.hh"
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
  
  auto eps = [](scalar_t r, scalar_t z) {
    return 1.78;
  };

  auto impulse_response = [](scalar_t t) {
    unsigned int order = 6;
    double tp = 2.0;   
    double retval = 1.0 / (tp * fact(order - 1)) * std::pow(t * (double)order / tp, (double)order) * std::exp(-t * (double)order / tp);    
    return retval;
  };
  
  CylinderGeometry geom(20, -20, 20, eps);
  InfEDipoleAntenna dipole(0.0, 10.0, -2.0, impulse_response);

  scalar_t t_end = 250;
  WeightingFieldCalculator wfc(geom, dipole, t_end);
  wfc.Calculate(10);
  
  return 0;
}
