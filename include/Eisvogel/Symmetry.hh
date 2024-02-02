#ifndef __SYMMETRY_HH
#define __SYMMETRY_HH

#include "CoordUtils.hh"

struct CylindricalSymmetry {

  static CoordVector GetOrbitIndex(const CoordVector& pos_txyz,
				   const CoordVector& start_coords_trz, const CoordVector& end_coords_trz,
				   const IndexVector& shape) {    
    CoordVector pos_trz = CoordUtils::TXYZ_to_TRZ(pos_txyz);
    return (pos_trz - start_coords_trz) / (end_coords_trz - start_coords_trz) * (CoordVector)shape;
  };
  
};

#endif
