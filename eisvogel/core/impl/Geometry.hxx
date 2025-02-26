#include "Vector.hh"

CylinderRegion::CylinderRegion() : CylinderRegion(0.0, 0.0, 0.0, 0.0) { }

CylinderRegion::CylinderRegion(scalar_t r_min, scalar_t r_max, scalar_t z_min, scalar_t z_max) :
  m_r_min(r_min), m_r_max(r_max), m_z_min(z_min), m_z_max(z_max) {
  throw_if_unphysical();
}

void CylinderRegion::SetRegion(scalar_t r_min, scalar_t r_max, scalar_t z_min, scalar_t z_max) {
  m_r_min = r_min;
  m_r_max = r_max;
  m_z_min = z_min;
  m_z_max = z_max;
  throw_if_unphysical();
}

void CylinderRegion::throw_if_unphysical() {
  if((m_r_min < (scalar_t)(0.0)) ||
     (m_r_max < (scalar_t)(0.0)) ||
     (m_r_min >= m_r_max)) {
    throw std::runtime_error("Error: unphysical configuration of CylinderGeometry");
  }
  if(m_z_min >= m_z_max) {
    throw std::runtime_error("Error: unphysical configuration of CylinderGeometry");
  }
}

CylinderGeometry::CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
				   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func) :
  m_r_max(r_max), m_z_min(z_min), m_z_max(z_max), m_eps_func(eps_func) {
  throw_if_unphysical();
}

void CylinderGeometry::throw_if_unphysical() {
  if((m_r_max < (scalar_t)(0.0)) || (m_z_min >= m_z_max)) {
    throw std::runtime_error("Error: unphysical configuration of CylinderGeometry");
  }
}

bool CylinderGeometry::contains(const CylinderRegion& reg) {
  if((reg.GetZMin() >= GetZMin()) &&
     (reg.GetZMax() < GetZMax()) &&
     (reg.GetRMin() >= GetRMin()) &&
     (reg.GetRMax() < GetRMax())) {
    return true;
  }
  return false;
}

double CylinderGeometry::eps(const meep::vec& pos) {
  RZCoordVector eisvogel_coords = toEisvogelCoords(pos);
  return m_eps_func(eisvogel_coords.r(), eisvogel_coords.z());
}

RZCoordVector CylinderGeometry::toEisvogelCoords(const meep::vec& meep_coords) {
  return RZCoordVector{(scalar_t)meep_coords.r(), (scalar_t)meep_coords.z() + m_z_min};
};

meep::vec CylinderGeometry::toMeepCoords(const RZCoordVector& coords) {
  return meep::veccyl(coords.r(), coords.z() - m_z_min);
}
