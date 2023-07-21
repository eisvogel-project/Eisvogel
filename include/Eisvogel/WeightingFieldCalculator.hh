#ifndef __WEIGHTING_FIELD_CALCULATOR__HH
#define __WEIGHTING_FIELD_CALCULATOR__HH

#include "Antenna.hh"
#include "Geometry.hh"
#include <meep.hpp>

class WeightingFieldCalculator {

public:
  WeightingFieldCalculator(const CylinderGeometry& geom, const Antenna& antenna,
			   double courant_factor = 0.5, double resolution = 20, double pml_width = 1.0);
  void Calculate();

private:
  std::shared_ptr<meep::grid_volume> gv;
  std::shared_ptr<meep::structure> s;
};

#endif
