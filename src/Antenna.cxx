#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

Antenna::Antenna(scalar_t z_pos, std::function<scalar_t(scalar_t)> impulse_response_func) :
  z_pos(z_pos), impulse_response_func(impulse_response_func) {
  is_integrated = false; // impulse_response_func directly specifies a current 
}

InfEDipoleAntenna::InfEDipoleAntenna(scalar_t z_pos, std::function<scalar_t(scalar_t)> impulse_response_func) :
  Antenna(z_pos, impulse_response_func) { };

void InfEDipoleAntenna::AddToGeometry(meep::fields& f, Geometry& geom) const {
  CoordVector antenna_pos = CoordUtils::MakeCoordVectorTRZ(0.0, 0.0, z_pos);
  
  // orientation hard-coded for now ... need to change
  meep::vec meeppos = geom.toMeepCoords(antenna_pos);
  std::cout << "meep pos" << std::endl;
  std::cout << "r = " << meeppos.r() << std::endl;
  std::cout << "z = " << meeppos.z() << std::endl;
  
  f.add_point_source(meep::Ez, *this, geom.toMeepCoords(antenna_pos));
}
