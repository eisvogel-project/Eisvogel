#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"
#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

int main(int argc, char* argv[]) {

  Geometry geom(20, -20, 20);
  Antenna dipole(-2.0);
  
  WeightingFieldCalculator wfc(geom, dipole);
  
  return 0;
}
