#include <iostream>
#include <fstream>
#include "Common.hh"
#include "WeightingFieldUtils.hh"
#include "Kernels.hh"
#include "Interpolator.hh"
#include "MathUtils.hh"

namespace WFU = WeightingFieldUtils;

int main(void) {

  CoordVector start_coords = CoordUtils::MakeCoordVectorTRZ(-0.2234, -2.0, -2.00876);
  CoordVector end_coords = CoordUtils::MakeCoordVectorTRZ(7.0, 3.4, 3.7);
  scalar_t tp = 5.0;
  unsigned int N = 4;
  scalar_t r_min = 0.1;

  std::string path = "/home/windischhofer/Eisvogel/Eisvogel/signal/build/electric_dipole_wf.bin";
  std::fstream ofs;
  
  ofs.open(path, std::ios::in | std::ios::binary);
  stor::Serializer ser(ofs);

  // WeightingField wf = ser.deserialize<WeightingField>();
  WeightingField wf = WFU::CreateElectricDipoleWeightingField(start_coords, end_coords, tp, N, r_min, 70);

  SplineInterpolationKernelOrder3 kernel;
  // SincInterpolationKernel kernel;
  Interpolator itpl_E_r(wf.E_r(), kernel);
  Interpolator itpl_E_z(wf.E_z(), kernel);

  CoordVector test_pos = CoordUtils::MakeCoordVectorTRZ(2.3, 0.27, 1.0);
  // CoordVector test_pos = CoordUtils::MakeCoordVectorTRZ(2.3, 0.27, 1.1);
  CoordVector test_frac_inds = wf.getFracInds(test_pos);
  
  // std::cout << "E_r = " << itpl_E_r.Interpolate(test_frac_inds) << std::endl;
  std::cout << "E_z = " << itpl_E_z.Interpolate(test_frac_inds) << std::endl;

  scalar_t Qw = 1.0;
  scalar_t eps0 = 1.0;  // vacuum dielectric constant
  scalar_t n = 1.0;  // refractive index of medium
  scalar_t ds = 1.0;
  scalar_t c = 1.0;  // speed of light in vacuum

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

  auto E_r = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
    r_xy = std::fabs(r_xy);
    scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
    scalar_t t_prop = r * n / c, t_del = t - t_prop;
    scalar_t cos_theta = z / r;
    
      return -2.0 * Qw * ds / (eps0 * 4 * M_PI) * cos_theta / std::pow(r, 3) * (filtered_theta(t_del, tp, N) + 
										t_prop * filtered_delta(t_del, tp, N));
  };

    auto E_theta = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t t_prop = r * n / c, t_del = t - t_prop;
      scalar_t sin_theta = r_xy / r;
      return -Qw * ds / (eps0 * 4 * M_PI) * sin_theta / std::pow(r, 3) * (filtered_theta(t_del, tp, N) + t_prop * filtered_delta(t_del, tp, N)
									  + std::pow(t_prop, 2) * filtered_delta_prime(t_del, tp, N));
    };

    auto E_z = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(std::pow(r_xy, 2) + std::pow(z, 2));
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * cos_theta - E_theta(t, r_xy, z) * sin_theta;
    };

    std::cout << "E_z(true) = " << E_z(CoordUtils::getT(test_pos), CoordUtils::getR(test_pos), CoordUtils::getZ(test_pos)) << std::endl;

  return 0;
}
