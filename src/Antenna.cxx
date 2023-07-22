#include "Eisvogel/Antenna.hh"
#include "Eisvogel/CoordUtils.hh"

unsigned int fact(unsigned arg) {
  unsigned int retval = 1;

  for(unsigned int cur = 1; cur <= arg; cur++) {
    retval *= cur;
  }
  
  return retval;
}

std::complex<double> srcfunc(double t, void*) {

  std::cout << "srcfunc hardcoded" << std::endl;
  
  unsigned int order = 6;
  double tp = 2.0;
  
  if(t > 0) {
    double retval = 1.0 / (tp * fact(order - 1)) * std::pow(t * (double)order / tp, (double)order) * std::exp(-t * (double)order / tp);
    std::cout << retval << std::endl;
    return retval;
  }
  return 0.0;
}

Antenna::Antenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
		 std::function<scalar_t(scalar_t)> impulse_response_func) :
  start_time(start_time), end_time(end_time), z_pos(z_pos), impulse_response_func(impulse_response_func) {
  is_integrated = false; // impulse_response_func directly specifies a current 
}

InfEDipoleAntenna::InfEDipoleAntenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
				     std::function<scalar_t(scalar_t)> impulse_response_func) :
  Antenna(start_time, end_time, z_pos, impulse_response_func) { };

void InfEDipoleAntenna::AddToGeometry(meep::fields& f, Geometry& geom) const {
  CoordVector antenna_pos = CoordUtils::MakeCoordVectorTRZ(0.0, 0.0, z_pos);
  
  // orientation hard-coded for now ... need to change
  meep::vec meeppos = geom.toMeepCoords(antenna_pos);
  std::cout << "meep pos" << std::endl;
  std::cout << "r = " << meeppos.r() << std::endl;
  std::cout << "z = " << meeppos.z() << std::endl;

  f.add_point_source(meep::Ez, *this, geom.toMeepCoords(antenna_pos));
  
  // meep::custom_src_time src(srcfunc, NULL);
  // f.add_point_source(meep::Ez, src, geom.toMeepCoords(antenna_pos));
}
