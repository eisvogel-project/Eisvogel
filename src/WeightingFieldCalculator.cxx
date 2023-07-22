#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  f = std::make_shared<meep::fields>(s.get());

  f -> output_hdf5(meep::Dielectric, gv -> surroundings());
  antenna.AddToGeometry(*f, geom);

  // =========
  
  while (f -> time() < 10) {
    std::cout << f -> time() << std::endl;
    f -> step();
  }

  f -> output_hdf5(meep::Ez, gv -> surroundings());
}
