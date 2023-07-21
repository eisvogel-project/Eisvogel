#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

unsigned int fact(unsigned arg) {
  unsigned int retval = 1;

  for(unsigned int cur = 1; cur <= arg; cur++) {
    retval *= cur;
  }
  
  return retval;
}

std::complex<double> srcfunc(double t, void*) {

  std::cout << "srcfunc hardcoded" << std::endl;
  
  unsigned int order = 6;
  double tp = 1.0;
  
  if(t > 0) {
    double retval = 1.0 / (tp * fact(order - 1)) * std::pow(t * (double)order / tp, (double)order) * std::exp(-t * (double)order / tp);
    return retval;
  }
  return 0.0;
}

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  f = std::make_shared<meep::fields>(s.get());

  f -> output_hdf5(meep::Dielectric, gv -> surroundings());

  meep::custom_src_time src(srcfunc, NULL);
  f -> add_point_source(meep::Ez, src, meep::veccyl(0.0, 15.0));

  // antenna.AddToGeometry(*f, geom);

  // =========
  
  while (f -> time() < 10) {
    std::cout << f -> time() << std::endl;
    f -> step();
  }

  f -> output_hdf5(meep::Ez, gv -> surroundings());
}
