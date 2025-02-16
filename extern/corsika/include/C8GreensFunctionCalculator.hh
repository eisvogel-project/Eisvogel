#pragma once

#include <functional>
#include <filesystem>

namespace GreensFunctionCalculator::MEEP {
  class CylindricalGreensFunctionCalculator;
}

class C8CylindricalGreensFunctionCalculator {

public:

  using eps_signature = float(float, float); // (r, z) -> eps
  using impulse_response_signature = float(float); // t -> current
  
  C8CylindricalGreensFunctionCalculator(float r_max, float z_min, float z_max, float z_ant,
					std::function<eps_signature> eps,
					std::function<impulse_response_signature> impulse_response);

  void calculate(std::filesystem::path gf_path, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir);
  
private:

  std::unique_ptr<GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator> m_gfc;

};
