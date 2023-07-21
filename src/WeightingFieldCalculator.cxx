#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(const CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {

  auto eps_meep = [geom](const meep::vec& pos) -> double {
    return geom.eps(pos.r(), pos.z());
  };
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(40, 40, resolution));
  s = std::make_shared<meep::structure>(*gv, eps_meep, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  
}
