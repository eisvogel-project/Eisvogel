#ifndef __COORD_UTILS_HH
#define __COORD_UTILS_HH

#include "Common.hh"
#include "NDArray.hh"

using CoordVector = DenseVector<scalar_t>;
using DeltaVector = DenseVector<scalar_t>;

namespace CoordUtils {

  static inline CoordVector MakeCoordVectorTRZ(scalar_t t, scalar_t r, scalar_t z) {
    return CoordVector({t, z, r});
  };
  static inline CoordVector MakeCoordVectorTXYZ(scalar_t t, scalar_t x, scalar_t y, scalar_t z) {
    return CoordVector({t, z, x, y});
  };
  static inline DeltaVector MakeDeltaVectorTRZ(scalar_t delta_t, scalar_t delta_r, scalar_t delta_z) {
    return DeltaVector({delta_t, delta_z, delta_r});
  };
  
  static inline scalar_t getT(const CoordVector& vect) {return vect(0);};
  static inline scalar_t getR(const CoordVector& vect) {return vect(2);};
  static inline scalar_t getX(const CoordVector& vect) {return vect(2);};
  static inline scalar_t getY(const CoordVector& vect) {return vect(3);};
  static inline scalar_t getZ(const CoordVector& vect) {return vect(1);};

  static inline std::size_t getTInd(const IndexVector& vect) {return vect(0);};
  static inline std::size_t getRInd(const IndexVector& vect) {return vect(2);};
  static inline std::size_t getXInd(const IndexVector& vect) {return vect(2);};
  static inline std::size_t getYInd(const IndexVector& vect) {return vect(3);};
  static inline std::size_t getZInd(const IndexVector& vect) {return vect(1);};
};

#endif
