#ifndef __ANTENNA_HH
#define __ANTENNA_HH

#include <functional>
#include <meep.hpp>
#include "CoordUtils.hh"
#include "Geometry.hh"

class Antenna : public meep::src_time {

public:
  Antenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
	  std::function<scalar_t(scalar_t)> impulse_response_func);
  virtual ~Antenna() { };
  
  virtual void AddToGeometry(meep::fields& f, Geometry& geom) const = 0;

public:
  virtual std::complex<double> current(double time, double dt) const {
    return dipole(time);
  };

  virtual std::complex<double> dipole(double time) const {
    if(time >= start_time && time <= end_time) {
      return impulse_response_func(time);
    } else {
      return 0.0;
    }
  };

  virtual double last_time() const {
    return end_time;
  };

protected:
  scalar_t start_time, end_time, z_pos;
  
private:
  std::function<scalar_t(scalar_t)> impulse_response_func;
  
};

class InfEDipoleAntenna : public Antenna {
  
public:
  InfEDipoleAntenna(scalar_t start_time, scalar_t end_time, scalar_t z_pos,
		    std::function<scalar_t(scalar_t)> impulse_response_func);
  virtual ~InfEDipoleAntenna() { };
  
  virtual void AddToGeometry(meep::fields& f, Geometry& geom) const;

  virtual InfEDipoleAntenna* clone() const { return new InfEDipoleAntenna(*this); }
  
};

#endif
