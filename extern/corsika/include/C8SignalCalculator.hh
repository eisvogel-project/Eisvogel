#pragma once

#include <filesystem>

class CylindricalGreensFunction;

class C8SignalCalculator {

public:
  C8SignalCalculator(std::filesystem::path gf_path);

  int calculate(float inval);
  
private:
  std::shared_ptr<CylindricalGreensFunction> m_gf;
  
};
