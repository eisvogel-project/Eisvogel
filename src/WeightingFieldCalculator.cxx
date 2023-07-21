#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(40, 40, resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  
}
