#include "WeightingFieldUtils.hh"

namespace WeightingFieldUtils{

  DenseWeightingField CreateElectricDipoleWeightingField() {
    
    std::size_t pts_t = 10, pts_r = 10, pts_z = 10;
    float t_min = 0.0, t_max = 1.0;
    float r_min = 0.0, r_max = 1.0;
    float z_min = 0.0, z_max = 1.0;
    
    auto E_r = [&](float t, float r, float z) -> float {
      return 0.0;
    };
    
    auto E_z = [&](float t, float r, float z) -> float {
      return 1.0;
    };
    
    auto E_phi = [&](float t, float r, float z) -> float {
      return 2.0;
    };
    
    DenseNDArray<float, 3> E_r_sampled({pts_t, pts_r, pts_z}, 0.0);
    DenseNDArray<float, 3> E_z_sampled({pts_t, pts_r, pts_z}, 0.0);
    DenseNDArray<float, 3> E_phi_sampled({pts_t, pts_r, pts_z}, 0.0);
    
    for(auto ind_t = 0; ind_t < pts_t; ind_t++) {
      for(auto ind_r = 0; ind_r < pts_r; ind_r++) {
	for(auto ind_z = 0; ind_z < pts_z; ind_z++) {

	  float t = (t_max - t_min) * ind_t / pts_t + t_min;
	  float r = (r_max - r_min) * ind_r / pts_r + r_min;
	  float z = (z_max - z_min) * ind_z / pts_z + z_min;

	  E_r_sampled(ind_t, ind_r, ind_z) = E_r(t, r, z);
	  E_z_sampled(ind_t, ind_r, ind_z) = E_z(t, r, z);
	  E_phi_sampled(ind_t, ind_r, ind_z) = E_phi(t, r, z);
	}
      }
    }
    
    return DenseWeightingField(std::move(E_r_sampled), std::move(E_z_sampled), std::move(E_phi_sampled),
			       t_min, t_max, r_min, r_max, z_min, z_max);
  }
}
