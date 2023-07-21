#ifndef __ANTENNA_HH
#define __ANTENNA_HH

#include <functional>
#include <meep.hpp>
#include "CoordUtils.hh"
#include "Geometry.hh"

class Antenna : public meep::src_time {

public:
  Antenna(scalar_t z_pos, std::function<scalar_t(scalar_t)> impulse_response_func);
  virtual void AddToGeometry(meep::fields& f, Geometry& geom) const = 0;

public:
  virtual std::complex<double> current(double time, double dt) const {
    return dipole(time);
  };

  virtual std::complex<double> dipole(double time) const {
    std::cout << "callback" << std::endl;
    
    if(time > 0.0) {
      return impulse_response_func(time);
    } else {
      return 0.0;
    }
  };

  virtual double last_time() const {
    std::cout << "last_time" << std::endl;
    return 10.0;
  };
  
private:
  std::function<scalar_t(scalar_t)> impulse_response_func;

protected:
    scalar_t z_pos;
  
};

class InfEDipoleAntenna : public Antenna {
  
public:
  InfEDipoleAntenna(scalar_t z_pos,  std::function<scalar_t(scalar_t)> impulse_response_func);
  virtual void AddToGeometry(meep::fields& f, Geometry& geom) const;
  
};

#endif
