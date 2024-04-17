#include "Eisvogel/GeometryOld.hh"

CylinderGeometry::CylinderGeometry(scalar_t r_max, scalar_t z_min, scalar_t z_max,
				   std::function<scalar_t(scalar_t r, scalar_t z)> eps_func) :
  r_max(r_max), z_min(z_min), z_max(z_max), eps_func(eps_func) { };

double CylinderGeometry::eps(const meep::vec& pos) {
  CoordVector eisvogel_coords = toEisvogelCoords(pos);
  return eps_func(CoordUtils::getR(eisvogel_coords), CoordUtils::getZ(eisvogel_coords));
}
