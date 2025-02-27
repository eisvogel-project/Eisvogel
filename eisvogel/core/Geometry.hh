#pragma once

#include <functional>
#include <meep.hpp>
#include "Common.hh"

struct RZCoordVector;

class CylinderRegion {

public:

  CylinderRegion();
  CylinderRegion(scalar_t r_min, scalar_t r_max, scalar_t z_min, scalar_t z_max);
  void SetRegion(scalar_t r_min, scalar_t r_max, scalar_t z_min, scalar_t z_max);

  bool IsEmpty();
  
  scalar_t GetRMin() const {return m_r_min;};
  scalar_t GetRMax() const {return m_r_max;};
  scalar_t GetZMin() const {return m_z_min;};
  scalar_t GetZMax() const {return m_z_max;};

  scalar_t GetRSize() const {return m_r_max - m_r_min;};
  scalar_t GetZSize() const {return m_z_max - m_z_min;};
  
private:
  void throw_if_unphysical();
  scalar_t m_r_min, m_r_max, m_z_min, m_z_max;
  
};

class CylinderGeometry : public meep::material_function {

public:

  CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
		   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func);

  scalar_t GetRMin() {return 0.0;};
  scalar_t GetRMax() {return m_r_max;};
  scalar_t GetZMin() {return m_z_min;};
  scalar_t GetZMax() {return m_z_max;};

  scalar_t GetRSize() {return m_r_max;};
  scalar_t GetZSize() {return m_z_max - m_z_min;};

  bool contains(const CylinderRegion& reg);
  
  meep::vec toMeepCoords(const RZCoordVector& coords);
  RZCoordVector toEisvogelCoords(const meep::vec& meep_coords);
  
public:
  
  virtual double eps(const meep::vec& pos);

  virtual double chi1p1(meep::field_type ft, const meep::vec &r) {
    (void)ft;
    return eps(r);
  }
  
private:
  void throw_if_unphysical();
  scalar_t m_r_max, m_z_min, m_z_max;

public:
  std::function<scalar_t(scalar_t r, scalar_t z)> m_eps_func;
  
};

#include "Geometry.hxx"
