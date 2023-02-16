#include <iostream>
#include "NDArray.hh"
#include "WeightingFieldUtils.hh"

namespace WFU = WeightingFieldUtils;

int main(void) {

  unsigned int len_t = 10;

  int n = 10;

  DenseWeightingField wf = WFU::CreateElectricDipoleWeightingField();

  std::cout << wf.E_r(0,0,0) << std::endl;
  std::cout << wf.E_z(0,0,0) << std::endl;
  std::cout << wf.E_phi(0,0,0) << std::endl;

  return 0;
}
