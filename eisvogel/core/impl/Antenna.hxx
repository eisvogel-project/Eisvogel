#include "Vector.hh"

Antenna::Antenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
		 std::function<scalar_t(scalar_t)> impulse_response_func) :
  m_start_time(start_time), m_end_time(end_time), m_z_pos(z_pos), m_scale_factor(1.0f),
  m_impulse_response_func(impulse_response_func) {

  // This is important: `is_integrated = false` specifies that MEEP expects us to
  // provide the feed current and not the dipole moment of the antenna
  is_integrated = false;
};

void Antenna::SetScaleFactor(scalar_t scale_factor) {
  m_scale_factor = scale_factor;
}

std::complex<double> Antenna::current(double time, [[maybe_unused]] double dt) const {
  
  float rtime = float(time);
  if(rtime >= m_start_time && rtime <= m_end_time) {
    return m_scale_factor * m_impulse_response_func(time) * (-1.0);  // Here is where the negative sign of the Green's function enters
  } else {
    return 0.0;
  }  
};

std::complex<double> Antenna::dipole([[maybe_unused]] double time) const {
  // Because we set `is_integrated = false`, no need to provide anything here.
  // Return NaN instead of 0.0 (or any other value) to make it obvious if this makes its way into any downstream calculations
  return std::numeric_limits<double>::quiet_NaN();
};

InfEDipoleAntenna::InfEDipoleAntenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
				     std::function<scalar_t(scalar_t)> impulse_response_func) :
  Antenna(start_time, end_time, z_pos, impulse_response_func) { };

void InfEDipoleAntenna::AddToGeometry(meep::fields& f, CylinderGeometry& geom) {

  // According to workaround described in https://github.com/NanoComp/meep/issues/2704, need to place source at r = 1.5 * delta_R,
  // where delta_R is the radial voxel extent
  scalar_t r_pos = 1.5 / f.a;
  RZCoordVector antenna_pos{r_pos, m_z_pos};
  
  // orientation hard-coded for now ... need to change
  f.add_point_source(meep::Ez, *this, geom.toMeepCoords(antenna_pos));
}
