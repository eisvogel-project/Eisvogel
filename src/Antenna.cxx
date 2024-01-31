#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

Antenna::Antenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
		 std::function<scalar_t(scalar_t)> impulse_response_func) :
  start_time(start_time), end_time(end_time), z_pos(z_pos), impulse_response_func(impulse_response_func) {

  // This is important: `is_integrated = false` specifies that MEEP expects us to
  // provide the feed current and not the dipole moment of the antenna
  is_integrated = false;
};

std::complex<double> Antenna::current(double time, double dt) const {
  
  float rtime = float(time);
  if(rtime >= start_time && rtime <= end_time) {
    return impulse_response_func(time);
  } else {
    return 0.0;
  }  
};

 std::complex<double> Antenna::dipole(double time) const {
   // Because we set `is_integrated = false`, no need to provide anything here.
   // Return NaN instead of 0.0 (or any other value) to make it obvious if this makes its way into any downstream calculations
   return std::numeric_limits<double>::quiet_NaN();
};

InfEDipoleAntenna::InfEDipoleAntenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
				     std::function<scalar_t(scalar_t)> impulse_response_func) :
  Antenna(start_time, end_time, z_pos, impulse_response_func) { };

void InfEDipoleAntenna::AddToGeometry(meep::fields& f, Geometry& geom) const {

  // According to workaround described in https://github.com/NanoComp/meep/issues/2704, need to place source at r = 1.5 * delta_R,
  // where delta_R is the radial voxel extent
  scalar_t r_pos = 1.5 / f.a;
  CoordVector antenna_pos = CoordUtils::MakeCoordVectorTRZ(0.0, r_pos, z_pos);
  
  // orientation hard-coded for now ... need to change
  f.add_point_source(meep::Ez, *this, geom.toMeepCoords(antenna_pos));
}
