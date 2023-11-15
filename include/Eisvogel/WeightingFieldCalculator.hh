#ifndef __WEIGHTING_FIELD_CALCULATOR__HH
#define __WEIGHTING_FIELD_CALCULATOR__HH

#include <memory>
#include <string>
#include <meep.hpp>
#include "Antenna.hh"
#include "Geometry.hh"

class WeightingFieldCalculator {

public:
  WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
			   double courant_factor = 0.5, double resolution = 20, double pml_width = 1.0);
  void Calculate(double t_end, std::string tmpdir = "");

private:
  std::shared_ptr<meep::grid_volume> gv;
  std::shared_ptr<meep::structure> s;
  std::shared_ptr<meep::fields> f;
};

#endif
