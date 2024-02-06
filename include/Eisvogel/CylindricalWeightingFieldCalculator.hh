#ifndef __CYLINDRICAL_WEIGHTING_FIELD_CALCULATOR__HH
#define __CYLINDRICAL_WEIGHTING_FIELD_CALCULATOR__HH

#include <memory>
#include <string>
#include <filesystem>
#include <meep.hpp>
#include "Antenna.hh"
#include "Geometry.hh"

class CylindricalWeightingFieldCalculator {

public:
  CylindricalWeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna, scalar_t t_end,
				      double courant_factor = 0.5, double resolution = 40, double pml_width = 1.0);
  void Calculate(std::filesystem::path outdir, std::filesystem::path tmpdir = "");

private:

  scalar_t m_t_end = 0.0;
  std::shared_ptr<CoordVector> m_start_coords;
  std::shared_ptr<CoordVector> m_end_coords;
  
  std::shared_ptr<meep::grid_volume> m_gv;
  std::shared_ptr<meep::structure> m_s;
  std::shared_ptr<meep::fields> m_f;
};

#endif
