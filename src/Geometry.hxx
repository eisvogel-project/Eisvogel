#include "Eisvogel/Geometry.hh"
#include "Vector.hh"

CylinderGeometry::CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
				   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func) :
  m_r_max(r_max), m_z_min(z_min), m_z_max(z_max), m_eps_func(eps_func) { }

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
