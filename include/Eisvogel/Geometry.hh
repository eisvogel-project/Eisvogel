#ifndef __GEOMETRY__HH
#define __GEOMETRY__HH

#include <functional>
#include <meep.hpp>
#include "Common.hh"
#include "CoordUtils.hh"

// class CartesianGeometry {

// };

// class CartesianAnisotropicGeometry {

// };

class CylinderGeometry : public meep::material_function {

public:

  CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
		   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func);
  
  scalar_t GetRMax() {return r_max;};
  scalar_t GetZMin() {return z_min;};
  scalar_t GetZMax() {return z_max;};

public:

  // Here are functions that meep calls
  virtual double eps(const meep::vec& pos);
  
private:
  scalar_t r_max, z_min, z_max;

public:
  std::function<scalar_t(scalar_t r, scalar_t z)> eps_func;
  
};

#endif
