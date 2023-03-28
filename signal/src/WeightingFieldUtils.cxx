#include "WeightingFieldUtils.hh"
#include "Common.hh"
#include "IteratorUtils.hh"
#include "CoordUtils.hh"
#include "MathUtils.hh"
#include <cmath>

namespace C = CoordUtils;

namespace WeightingFieldUtils {

  WeightingField CreateElectricDipoleWeightingField(const CoordVector& start_coords, const CoordVector& end_coords,
						    scalar_t tp, unsigned int N, scalar_t os_factor) {

    scalar_t Qw = 1.0;
    scalar_t eps0 = 1.0;  // vacuum dielectric constant
    scalar_t n = 1.0;  // refractive index of medium
    scalar_t ds = 1.0;
    scalar_t c = 1.0;  // speed of light in vacuum

    // compute required step size for sampling of weighting field
    scalar_t fmax = (scalar_t)N / (2 * M_PI * tp) * std::sqrt(std::pow(2.0, 1.0 / (N + 1)) - 1);
    scalar_t lambda_min = c / (fmax * n);
    scalar_t delta_t = 1.0 / (2 * fmax * os_factor);
    scalar_t delta_pos = lambda_min / (2.0 * os_factor);
    
    std::cout << "---------------------------" << std::endl;
    std::cout << "Using oversampling factor = " << os_factor << std::endl;
    std::cout << "delta_t = " << delta_t << std::endl;
    std::cout << "delta_r = delta_z = " << delta_pos << std::endl;
    std::cout << "---------------------------" << std::endl;

    DeltaVector step = C::MakeCoordVectorTRZ(delta_t, delta_pos, delta_pos);

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

    // Weighting field in spherical coordinates
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

    // Weighting field in cylindrical coordinates
    auto E_rxy = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * sin_theta + E_theta(t, r_xy, z) * cos_theta;
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

      scalar_t cur_E_rxy = E_rxy(t, r, z);
      scalar_t cur_E_z = E_z(t, r, z);
      scalar_t cur_E_phi = E_phi(t, r, z);
      
      if(!(std::isfinite(cur_E_rxy) && std::isfinite(cur_E_z) && std::isfinite(cur_E_phi))) {
	throw;
      }

      E_r_sampled(ind) = cur_E_rxy;
      E_z_sampled(ind) = cur_E_z;
      E_phi_sampled(ind) = cur_E_phi;
    }
    
    return WeightingField(std::move(E_r_sampled), std::move(E_z_sampled), std::move(E_phi_sampled),
			  std::move(start_coords), std::move(end_coords));
  }
}
