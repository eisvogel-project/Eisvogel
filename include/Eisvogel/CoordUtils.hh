#ifndef __COORD_UTILS_HH
#define __COORD_UTILS_HH

#include "Common.hh"
#include "NDArray.hh"
#include <cmath>

using CoordVector = DenseVector<scalar_t>;
using DeltaVector = DenseVector<scalar_t>;
using FieldVector = DenseVector<scalar_t>;

namespace CoordUtils {

  static inline CoordVector MakeCoordVectorTRZ(scalar_t t, scalar_t r, scalar_t z) {
    return CoordVector({t, z, r});
  }
  static inline CoordVector MakeCoordVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z) {
    return CoordVector({t, z, x, y});
  }
  static inline DeltaVector MakeDeltaVectorTRZ(scalar_t delta_t, scalar_t delta_r, scalar_t delta_z) {
    return DeltaVector({delta_t, delta_z, delta_r});
  }
  static inline DeltaVector MakeDeltaVectorTXYZ(scalar_t delta_t, scalar_t delta_x, scalar_t delta_y, scalar_t delta_z) {
    return DeltaVector({delta_t, delta_z, delta_x, delta_y});
  }
  
  static inline scalar_t getT(const CoordVector& vect) {return vect(0);}
  static inline scalar_t getR(const CoordVector& vect) {return vect(2);}
  static inline scalar_t getX(const CoordVector& vect) {return vect(2);}
  static inline scalar_t getY(const CoordVector& vect) {return vect(3);}
  static inline scalar_t getZ(const CoordVector& vect) {return vect(1);}

  static inline std::size_t getTInd(const IndexVector& vect) {return vect(0);}
  static inline std::size_t getRInd(const IndexVector& vect) {return vect(2);}
  static inline std::size_t getXInd(const IndexVector& vect) {return vect(2);}
  static inline std::size_t getYInd(const IndexVector& vect) {return vect(3);}
  static inline std::size_t getZInd(const IndexVector& vect) {return vect(1);}

  static inline CoordVector TXYZ_to_TRZ(const CoordVector& point_txyz) {
    return MakeCoordVectorTRZ(getT(point_txyz),
			      std::sqrt(std::pow(getX(point_txyz), 2) + std::pow(getY(point_txyz), 2)),
			      getZ(point_txyz)
			      );
  }

  static inline FieldVector MakeFieldVectorRZPHI(scalar_t comp_r, scalar_t comp_z, scalar_t comp_phi) {
    return CoordVector({comp_r, comp_z, comp_phi});
  }
  static inline FieldVector MakeFieldVectorXYZ(scalar_t comp_x, scalar_t comp_y, scalar_t comp_z) {
    return CoordVector({comp_x, comp_z, comp_y});
  }
  static inline scalar_t getXComponent(const FieldVector& vect) {return vect(0);}
  static inline scalar_t getYComponent(const FieldVector& vect) {return vect(2);}
  static inline scalar_t getZComponent(const FieldVector& vect) {return vect(1);}
  static inline scalar_t getRComponent(const FieldVector& vect) {return vect(0);}
  static inline scalar_t getPHIComponent(const FieldVector& vect) {return vect(2);}

  static inline FieldVector RZPHI_to_XYZ(const FieldVector& field_rzphi, const CoordVector& point_txyz) {
    scalar_t r_xy = std::sqrt(std::pow(getX(point_txyz), 2) + std::pow(getY(point_txyz), 2));
    scalar_t cos_phi = 0.0, sin_phi = 0.0;
    if(r_xy > 0) {
      cos_phi = getX(point_txyz) / r_xy;
      sin_phi = getY(point_txyz) / r_xy;
    }

    return MakeFieldVectorXYZ(cos_phi * getRComponent(field_rzphi) - sin_phi * getPHIComponent(field_rzphi),
			      sin_phi * getRComponent(field_rzphi) + cos_phi * getPHIComponent(field_rzphi),
			      getZComponent(field_rzphi));
  }

};

#endif
