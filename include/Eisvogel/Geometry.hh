#ifndef __GEOMETRY__HH
#define __GEOMETRY__HH

#include <functional>
#include <meep.hpp>
#include "Common.hh"
#include "CoordUtils.hh"

class Geometry : public meep::material_function {

public:
  virtual double eps(const meep::vec& pos) = 0;
  virtual meep::vec toMeepCoords(const CoordVector& coords) = 0;
  virtual double chi1p1(meep::field_type ft, const meep::vec &r) {
    (void)ft;
    return eps(r);
  }
};

class CylinderGeometry : public Geometry {

public:

  CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
		   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func);
  
  scalar_t GetRMax() {return r_max;};
  scalar_t GetZMin() {return z_min;};
  scalar_t GetZMax() {return z_max;};

  meep::vec toMeepCoords(const CoordVector& coords) {
    return meep::veccyl(CoordUtils::getR(coords), CoordUtils::getZ(coords) - z_min);
  };
  
public:
  virtual double eps(const meep::vec& pos);
  
private:
  scalar_t r_max, z_min, z_max;

public:
  std::function<scalar_t(scalar_t r, scalar_t z)> eps_func;
  
};

// class CartesianGeometry {

// };

// class CartesianAnisotropicGeometry {

// };

#endif
