#ifndef __ANTENNA_HH
#define __ANTENNA_HH

#include "CoordUtils.hh"

class Antenna {

public: 
  Antenna(scalar_t z_pos);  
  void GetSource();

private:
  scalar_t z_pos; 
};

#endif
