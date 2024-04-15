#include "Eisvogel/GreensFunctionUtils.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "DistributedNDVecArray.hh"
#include "GreensFunction.hh"

namespace GreensFunctionUtils {

  void CreateElectricDipoleGreensFunction(std::filesystem::path gf_path, const RZCoordVector& start_coords, const RZCoordVector& end_coords, scalar_t t_end, scalar_t ior,
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

    scalar_t t_start = 0.0;
    
    std::cout << "---------------------------" << std::endl;
    std::cout << "Using oversampling factor = " << os_factor << std::endl;
    std::cout << "delta_t = " << delta_t << std::endl;
    std::cout << "delta_r = delta_z = " << delta_pos << std::endl;
    std::cout << "start_coords: t = " << t_start << ", r = " << start_coords.r() << ", z = " << start_coords.z() << std::endl;
    std::cout << "end_coords: t = " << t_end << ", r = " << end_coords.r() << ", z = " << end_coords.z() << std::endl;
    std::cout << "---------------------------" << std::endl;

    RZTVector<scalar_t> stepsize_requested{delta_pos, delta_pos, delta_t};
    RZTCoordVector start_coords_rzt;
    RZTCoordVector end_coords_rzt;
    
    auto fill_and_register_chunk = [](const RZTIndexVector& chunk_start, const RZTIndexVector& chunk_end){

      
    };

    RZTIndexVector start_inds(0);
    RZTIndexVector end_inds = ((end_coords_rzt - start_coords_rzt) / stepsize_requested).template as_type<std::size_t>();
    RZTVector<std::size_t> chunk_size(max_pts_in_chunk);
    IteratorUtils::index_loop_over_chunks(start_inds, end_inds, chunk_size);
    
    // DeltaVector domain_size = end_coords - start_coords;
    // CoordVector number_pts = domain_size / stepsize_requested;

    // IndexVector number_chunks = number_pts / max_pts_in_chunk + 1;
    // DeltaVector chunk_size = domain_size / (DeltaVector)(number_chunks);
    // IndexVector number_pts_in_chunk = chunk_size / stepsize_requested + 1;    

  }
}
