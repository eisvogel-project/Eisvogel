#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

Antenna::Antenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
		 std::function<scalar_t(scalar_t)> impulse_response_func) :
  start_time(start_time), end_time(end_time), z_pos(z_pos), impulse_response_func(impulse_response_func) { };

std::complex<double> Antenna::current(double time, double dt) const {
  if (is_integrated)
    return src_time::current(time, dt);
  else
    return dipole(time);
};

std::complex<double> Antenna::dipole(double time) const {
  float rtime = float(time);
  if(rtime >= start_time && rtime <= end_time) {
    return impulse_response_func(time);
  } else {
    return 0.0;
  }
};

InfEDipoleAntenna::InfEDipoleAntenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
				     std::function<scalar_t(scalar_t)> impulse_response_func) :
  Antenna(start_time, end_time, z_pos, impulse_response_func) { };

void InfEDipoleAntenna::AddToGeometry(meep::fields& f, Geometry& geom) const {
  CoordVector antenna_pos = CoordUtils::MakeCoordVectorTRZ(0.0, 1, z_pos);
  
  // orientation hard-coded for now ... need to change
  f.add_point_source(meep::Ez, *this, geom.toMeepCoords(antenna_pos));
}
