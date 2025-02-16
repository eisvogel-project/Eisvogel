#pragma once

#include <functional>
#include <filesystem>

// Forward declarations
namespace GreensFunctionCalculator::MEEP {
  class CylindricalGreensFunctionCalculator;
}
class CylinderGeometry;
class Antenna;
namespace meep {
  class initialize;
}

class C8CylindricalGreensFunctionCalculator {

public:
  
  C8CylindricalGreensFunctionCalculator(int argc, char* argv[]);
  ~C8CylindricalGreensFunctionCalculator();

  using eps_signature = float(float, float); // (r, z) -> eps
  using impulse_response_signature = float(float); // t -> current
  
  void calculate(float r_max, float z_min, float z_max, float z_ant, float t_end,
		 std::function<eps_signature> eps,
		 std::function<impulse_response_signature> impulse_response,
		 std::filesystem::path gf_path, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir);
  
private:

  std::unique_ptr<meep::initialize> m_mpi;
  std::unique_ptr<GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator> m_gfc;
  std::unique_ptr<CylinderGeometry> m_geom;
  std::unique_ptr<Antenna> m_ant;

};
