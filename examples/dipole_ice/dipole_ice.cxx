#include <iostream>
#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingFieldCalculator.hh"
#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

int main(int argc, char* argv[]) {

  auto eps = [](scalar_t r, scalar_t z) {
    if (z < 0)
      return 1.78;
    return 1.0;
  };
  
  CylinderGeometry geom(20, -20, 20, eps);
  Antenna dipole(-2.0);
  
  WeightingFieldCalculator wfc(geom, dipole);
  
  return 0;
}
