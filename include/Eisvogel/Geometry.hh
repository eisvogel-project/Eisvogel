#pragma once

#include <functional>
#include <meep.hpp>
#include "Common.hh"

struct RZCoordVector;

class CylinderGeometry : public meep::material_function {

public:

  CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
		   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func);
  
  scalar_t GetRMax() {return m_r_max;};
  scalar_t GetZMin() {return m_z_min;};
  scalar_t GetZMax() {return m_z_max;};

  meep::vec toMeepCoords(const RZCoordVector& coords);
  RZCoordVector toEisvogelCoords(const meep::vec& meep_coords);
  
public:
  
  virtual double eps(const meep::vec& pos);

  virtual double chi1p1(meep::field_type ft, const meep::vec &r) {
    (void)ft;
    return eps(r);
  }
  
private:
  scalar_t m_r_max, m_z_min, m_z_max;

public:
  std::function<scalar_t(scalar_t r, scalar_t z)> m_eps_func;
  
};

// class CartesianGeometry {

// };

// class CartesianAnisotropicGeometry {

// };
