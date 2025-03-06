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
		 std::filesystem::path gf_path, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
		 double courant_factor = 0.5, double resolution = 28, double timestep = 0.1, double pml_width = 1.0, std::size_t downsampling_on_disk = 2,
		 double dynamic_range = 100, double abs_min_field = 1e-20, std::size_t chunk_overlap = 2, std::size_t chunk_size_linear = 800,
		 std::size_t rechunk_cache_depth = 5);
  
private:

  std::unique_ptr<meep::initialize> m_mpi;
  std::unique_ptr<GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator> m_gfc;
  std::unique_ptr<CylinderGeometry> m_geom;
  std::unique_ptr<Antenna> m_ant;

};
