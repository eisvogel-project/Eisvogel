#pragma once

#include <memory>
#include <filesystem>
#include <meep.hpp>
#include "Geometry.hh"
#include "Antenna.hh"

struct RZTCoordVector;

class CylindricalGreensFunctionCalculator {

public:

  CylindricalGreensFunctionCalculator(CylinderGeometry& geom, const Antenna& antenna, scalar_t t_end,
				      double courant_factor = 0.5, double resolution = 12, double pml_width = 1.0);

  void Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir);
  
private:

  scalar_t m_t_end;
  std::shared_ptr<RZTCoordVector> m_start_coords;
  std::shared_ptr<RZTCoordVector> m_end_coords;
  
  std::shared_ptr<meep::grid_volume> m_gv;
  std::shared_ptr<meep::structure> m_s;
  std::shared_ptr<meep::fields> m_f;  
};
