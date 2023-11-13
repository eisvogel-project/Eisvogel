#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  f = std::make_shared<meep::fields>(s.get());

  antenna.AddToGeometry(*f, geom);
}

void WeightingFieldCalculator::Calculate(double t_end) {

  f -> output_hdf5(meep::Dielectric, gv -> surroundings());
    
  while (f -> time() < t_end) {
    if(meep::am_master()) {
      std::cout << "Simulation time: " << f -> time() << std::endl;
    }
    f -> step();
    f -> output_hdf5(meep::Ez, gv -> surroundings());
  }
}
