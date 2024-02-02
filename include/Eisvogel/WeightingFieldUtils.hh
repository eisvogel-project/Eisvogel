#ifndef __WEIGHTING_FIELD_UTILS_HH
#define __WEIGHTING_FIELD_UTILS_HH

#include "CoordUtils.hh"
#include "NDArray.hh"

namespace WeightingFieldUtils {

  void CreateElectricDipoleWeightingField(std::string wf_path,
					  const CoordVector& start_coords, const CoordVector& end_coords,
					  scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor = 1.5, scalar_t n = 1,
					  unsigned int max_pts_in_chunk = 400);

  void SampleElectricDipoleWeightingFieldChunk(ScalarField3D<scalar_t>& E_r_buffer, ScalarField3D<scalar_t>& E_z_buffer, ScalarField3D<scalar_t>& E_phi_buffer,
					       const CoordVector& start_coords, const CoordVector& end_coords, 
					       scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor, scalar_t n);

}

#endif
