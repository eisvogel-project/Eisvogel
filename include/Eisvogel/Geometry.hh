#ifndef __GEOMETRY__HH
#define __GEOMETRY__HH

#include "Common.hh"

class Geometry {

public:
  Geometry(scalar_t r_max, scalar_t z_min, scalar_t z_max);

private:
  scalar_t r_max, z_min, z_max;  
};

#endif
