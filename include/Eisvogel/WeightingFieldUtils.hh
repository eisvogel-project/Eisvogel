#ifndef __WEIGHTING_FIELD_UTILS_HH
#define __WEIGHTING_FIELD_UTILS_HH

#include "WeightingField.hh"
#include "NDArray.hh"

namespace WeightingFieldUtils {

<<<<<<< HEAD
WeightingField CreateElectricDipoleWeightingField(const CoordVector& start_coords, const CoordVector& end_coords,
						  scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor = 1.5,
                                                  scalar_t refractive_index=1);
=======
  void CreateElectricDipoleWeightingField(std::string wf_path,
					  const CoordVector& start_coords, const CoordVector& end_coords,
					  scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor = 1.5);

  WeightingField SampleElectricDipoleWeightingField(const CoordVector& start_coords, const CoordVector& end_coords,
						    scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor);
>>>>>>> bde6b98337f452870b0f96be9d481e6a33acf76e

}

#endif
