#include "C8SignalCalculator.hh"

#include <iostream>

C8SignalCalculator::C8SignalCalculator(std::filesystem::path gf_path) {

  std::cout << "in interface constructor" << std::endl;  
}

int C8SignalCalculator::calculate(float inval) {
  return 18;
}
