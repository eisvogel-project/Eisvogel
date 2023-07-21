#ifndef __WEIGHTING_FIELD_CALCULATOR__HH
#define __WEIGHTING_FIELD_CALCULATOR__HH

#include "Antenna.hh"
#include "Geometry.hh"

class WeightingFieldCalculator {

public:
  WeightingFieldCalculator(const Geometry& geom, const Antenna& antenna,
			   double courant_factor = 0.5, double resolution = 20);
  void Calculate();

private:
  double resolution, courant_factor;
};

#endif
