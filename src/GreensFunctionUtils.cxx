#include "Eisvogel/GreensFunctionUtils.hh"
#include "DistributedNDVecArray.hh"
#include "GreensFunction.hh"

namespace GreensFunctionUtils {

  void CreateElectricDipoleGreensFunction(std::filesystem::path gf_path, const RZTCoordVector& start_coords, const RZTCoordVector& end_coords, scalar_t ior,
					  scalar_t filter_t_peak, unsigned int filter_order, scalar_t r_min,
					  scalar_t os_factor, std::size_t max_pts_in_chunk) {
    
    // make sure to start from scratch
    std::filesystem::remove_all(gf_path);
    
    // compute required step size for sampling of weighting field
    scalar_t c = 1.0;  // speed of light in vacuum
    scalar_t fmax = (scalar_t)(filter_order) / (2 * M_PI * filter_t_peak) * std::sqrt(std::pow(2.0, 1.0 / (filter_order + 1)) - 1);
    scalar_t lambda_min = c / (fmax * ior);
    scalar_t delta_t = 1.0 / (2 * fmax * os_factor);
    scalar_t delta_pos = lambda_min / (2.0 * os_factor);
    
    std::cout << "---------------------------" << std::endl;
    std::cout << "Using oversampling factor = " << os_factor << std::endl;
    std::cout << "delta_t = " << delta_t << std::endl;
    std::cout << "delta_r = delta_z = " << delta_pos << std::endl;
    std::cout << "start_coords: t = " << start_coords.t() << ", r = " << start_coords.r() << ", z = " << start_coords.z() << std::endl;
    std::cout << "---------------------------" << std::endl;    
  }
}
