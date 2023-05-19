#ifndef __WEIGHTING_FIELD_UTILS_HH
#define __WEIGHTING_FIELD_UTILS_HH

#include "WeightingField.hh"
#include "NDArray.hh"

namespace WeightingFieldUtils {

  void CreateElectricDipoleWeightingField(std::string wf_path,
					  const CoordVector& start_coords, const CoordVector& end_coords,
					  scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor = 1.5);

  WeightingField SampleElectricDipoleWeightingField(const CoordVector& start_coords, const CoordVector& end_coords,
						    scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor);

}

#endif
