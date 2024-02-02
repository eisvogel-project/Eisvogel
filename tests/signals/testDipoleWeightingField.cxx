#include <iostream>

#include "Eisvogel/WeightingField.hh"
#include "Eisvogel/Symmetry.hh"

int main(int argc, char* argv[]) {

  std::string wf_path = argv[1];
  
  WeightingField<TRZFieldIndexer<CylindricalSymmetry>, RZFieldStorage> bla(wf_path);

  bla.Er({1.0, 1.0, 1.0, 1.0});
  
  std::cout << "bla" << std::endl;
}
