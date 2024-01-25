#include "Eisvogel/WeightingFieldUtils.hh"
#include "Eisvogel/Common.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "Eisvogel/CoordUtils.hh"
#include "Eisvogel/MathUtils.hh"
#include "Eisvogel/Serialization.hh"
#include "Eisvogel/DistributedWeightingField.hh"
#include <cmath>
#include <iostream>
#include <fstream>

namespace C = CoordUtils;

namespace WeightingFieldUtils {

  void CreateElectricDipoleWeightingField(std::string wf_path,
					  const CoordVector& start_coords, const CoordVector& end_coords,
					  scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor, scalar_t n) {
    
    DistributedWeightingField dwf(wf_path, start_coords, end_coords);

    // compute required step size for sampling of weighting field
    scalar_t c = 1.0;  // speed of light in vacuum
    scalar_t fmax = (scalar_t)N / (2 * M_PI * tp) * std::sqrt(std::pow(2.0, 1.0 / (N + 1)) - 1);
    scalar_t lambda_min = c / (fmax * n);
    scalar_t delta_t = 1.0 / (2 * fmax * os_factor);
    scalar_t delta_pos = lambda_min / (2.0 * os_factor);
    
    std::cout << "---------------------------" << std::endl;
    std::cout << "Using oversampling factor = " << os_factor << std::endl;
    std::cout << "delta_t = " << delta_t << std::endl;
    std::cout << "delta_r = delta_z = " << delta_pos << std::endl;
    std::cout << "start_coords: t = " << C::getT(start_coords) << ", r = " << C::getR(start_coords) << ", z = " << C::getZ(start_coords) << std::endl;
    std::cout << "---------------------------" << std::endl;

    DeltaVector stepsize_requested = C::MakeCoordVectorTRZ(delta_t, delta_pos, delta_pos);
    
    CoordVector number_pts = (end_coords - start_coords) / stepsize_requested;
    std::size_t pts_t = C::getT(number_pts), pts_r = C::getR(number_pts), pts_z = C::getZ(number_pts);

    // TODO: generalize this to produce more than one chunk
    
    ScalarField3D<scalar_t> chunk_buffer_E_r({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> chunk_buffer_E_z({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> chunk_buffer_E_phi({pts_t, pts_z, pts_r}, 0.0);
    
    // Fill chunk buffers
    SampleElectricDipoleWeightingFieldChunk(chunk_buffer_E_r, chunk_buffer_E_z, chunk_buffer_E_phi, start_coords, end_coords, tp, N, r_min, os_factor, n);

    // Register chunk buffers
    dwf.RegisterChunk(chunk_buffer_E_r, chunk_buffer_E_z, chunk_buffer_E_phi, IndexVector({0, 0, 0}));

    dwf.Flush();
  }
  
  // TODO: three separate `ScalarField3D`s to be replaced with single vector field
  void SampleElectricDipoleWeightingFieldChunk(ScalarField3D<scalar_t>& E_r_buffer, ScalarField3D<scalar_t>& E_z_buffer, ScalarField3D<scalar_t>& E_phi_buffer,
					       const CoordVector& start_coords, const CoordVector& end_coords, 
					       scalar_t tp, unsigned int N, scalar_t r_min, scalar_t os_factor, scalar_t n) {

    scalar_t Qw = 1.0;
    scalar_t eps0 = 1.0;  // vacuum dielectric constant
    scalar_t ds = 1.0;
    scalar_t c = 1.0;  // speed of light in vacuum
    
    auto filtered_theta = [&](scalar_t t) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return 1.0 - MathUtils::incomplete_gamma(1 + N, N * t / tp) / std::exp(std::lgamma(N + 1));
    };
    
    auto filtered_delta = [&](scalar_t t) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return std::pow(t / tp * N, N) * std::exp(-t / tp * N) / (tp * std::exp(std::lgamma(N)));
    };

    auto filtered_delta_prime = [&](scalar_t t) -> scalar_t {
      if(t <= 0) {
	return 0.0;
      }
      return filtered_delta(t) * (tp - t) * N / (tp * t);
    };  

    // Weighting field in spherical coordinates
    auto E_r = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      if(r < r_min) {
	return std::nan("");
      }
      scalar_t t_prop = r * n / c, t_del = t - t_prop;
      scalar_t cos_theta = z / r;

      return -2.0 * Qw * ds / (eps0 * 4 * M_PI) * cos_theta / std::pow(r, 3) * (filtered_theta(t_del) + 
										t_prop * filtered_delta(t_del));
    };

    auto E_theta = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      if(r < r_min) {
	return std::nan("");
      }
      scalar_t t_prop = r * n / c, t_del = t - t_prop;
      scalar_t sin_theta = r_xy / r;
      return -Qw * ds / (eps0 * 4 * M_PI) * sin_theta / std::pow(r, 3) * (filtered_theta(t_del) + t_prop * filtered_delta(t_del)
									  + std::pow(t_prop, 2) * filtered_delta_prime(t_del));
    };

    // Weighting field in cylindrical coordinates
    auto E_rxy = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * sin_theta + E_theta(t, r_xy, z) * cos_theta;
    };
    
    auto E_z = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * cos_theta - E_theta(t, r_xy, z) * sin_theta;
    };
    
    auto E_phi = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      return 0.0;
    };

    IndexVector start_inds({0, 0, 0});
    IndexVector end_inds(E_r_buffer.shape());

    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      IndexVector ind = cnt.index();

      CoordVector coords = WeightingField::FracIndsToCoord(ind, start_coords, end_coords, E_r_buffer.shape());
      scalar_t t = C::getT(coords);
      scalar_t r = C::getR(coords);
      scalar_t z = C::getZ(coords);

      E_r_buffer(ind) = E_rxy(t, r, z);
      E_z_buffer(ind) = E_z(t, r, z);
      E_phi_buffer(ind) = E_phi(t, r, z);
      
    }
  }
}
