#include "WeightingFieldUtils.hh"
#include "Common.hh"
#include "IteratorUtils.hh"
#include "CoordUtils.hh"
#include "NDArray.hh"
#include "MathUtils.hh"
#include <cmath>

namespace C = CoordUtils;

namespace WeightingFieldUtils {

  WeightingField CreateElectricDipoleWeightingField() {
    
    // These will be arguments eventually
    CoordVector start_coords = C::MakeCoordVectorTRZ(0.0, 0.1, -10.0);
    CoordVector end_coords = C::MakeCoordVectorTRZ(320.0, 300.0, 30.0);
    DeltaVector step = C::MakeCoordVectorTRZ(0.1, 0.1, 1); // step size

    scalar_t Qw = 1.0;
    scalar_t eps0 = 1.0;
    scalar_t n = 1.0;
    scalar_t ds = 1.0;

    scalar_t tp = 5.0;
    unsigned int N = 1;

    // These will be stored in some central place
    scalar_t c = 1.0;

    auto filtered_theta = [&](scalar_t t, scalar_t tp, unsigned int N) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return 1.0 - MathUtils::incomplete_gamma(1 + N, N * t / tp) / std::exp(std::lgamma(N + 1));
    };

    auto filtered_delta = [&](scalar_t t, scalar_t tp, unsigned int N) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return std::pow(t / tp * N, N) * std::exp(-t / tp * N) / (tp * std::exp(std::lgamma(N)));
    };

    auto filtered_delta_prime = [&](scalar_t t, scalar_t tp, unsigned int N) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return filtered_delta(t, tp, N) * (tp - t) * N / (tp * t);
    };  
    
    // ==============

    CoordVector number_pts = (end_coords - start_coords) / step;
    std::size_t pts_t = C::getT(number_pts), pts_r = C::getR(number_pts), pts_z = C::getZ(number_pts);

    auto E_r = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t t_prop = r * n / c, t_del = t - t_prop;
      scalar_t cos_theta = z / r;

      return -2.0 * Qw * ds / (eps0 * 4 * M_PI) * cos_theta / std::pow(r, 3) * (filtered_theta(t_del, tp, N) + 
										t_prop * filtered_delta(t_del, tp, N));
    };

    auto E_theta = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t t_prop = r * n / c, t_del = t - t_prop;
      scalar_t sin_theta = r_xy / r;
      return -Qw * ds / (eps0 * 4 * M_PI) * sin_theta / std::pow(r, 3) * (filtered_theta(t_del, tp, N) + t_prop * filtered_delta(t_del, tp, N)
									  + std::pow(t_prop, 2) * filtered_delta_prime(t_del, tp, N));
    };
    
    auto E_z = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * cos_theta - E_theta(t, r_xy, z) * sin_theta;
    };
    
    auto E_phi = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      return 0.0;
    };
    
    // TODO: find a better way to make sure the ordering t, z, r etc is not messed up
    ScalarField3D<scalar_t> E_r_sampled({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> E_z_sampled({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> E_phi_sampled({pts_t, pts_z, pts_r}, 0.0);
    
    IndexVector start_inds({0, 0, 0});
    IndexVector end_inds({pts_t, pts_z, pts_r});

    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      IndexVector ind = cnt.index();

      scalar_t t = C::getT(start_coords) + C::getTInd(ind) * C::getT(step);
      scalar_t r = C::getR(start_coords) + C::getRInd(ind) * C::getR(step);
      scalar_t z = C::getZ(start_coords) + C::getZInd(ind) * C::getZ(step);

      scalar_t cur_E_r = E_r(t, r, z);
      scalar_t cur_E_z = E_z(t, r, z);
      scalar_t cur_E_phi = E_phi(t, r, z);
      
      if(!(std::isfinite(cur_E_r) && std::isfinite(cur_E_z) && std::isfinite(cur_E_phi))) {
	throw;
      }

      E_r_sampled(ind) = cur_E_r;
      E_z_sampled(ind) = cur_E_z;
      E_phi_sampled(ind) = cur_E_phi;
    }
    
    return WeightingField(std::move(E_r_sampled), std::move(E_z_sampled), std::move(E_phi_sampled),
			  std::move(start_coords), std::move(end_coords));
  }
}
