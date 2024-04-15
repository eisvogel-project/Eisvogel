#include <cassert>
#include "Eisvogel/GreensFunctionUtils.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "DistributedNDVecArray.hh"
#include "GreensFunction.hh"
#include "Symmetry.hh"

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

    RZTCoordVector start_coords_rzt{start_coords.r(), start_coords.z(), t_start};
    RZTCoordVector end_coords_rzt{end_coords.r(), end_coords.z(), t_end};

    // Determine step size so that an integer number of samples fits into the domain of the Green's function
    RZTVector<scalar_t> stepsize_requested{delta_pos, delta_pos, delta_t};
    RZTVector<std::size_t> number_pts = ((end_coords_rzt - start_coords_rzt) / stepsize_requested).template as_type<std::size_t>() + 1;
    RZTVector<scalar_t> stepsize = (end_coords_rzt - start_coords_rzt) / number_pts.template as_type<scalar_t>();
    
    // Resulting start- and end indices
    RZTIndexVector start_inds(0);
    RZTIndexVector end_inds = ((end_coords_rzt - start_coords_rzt) / stepsize).template as_type<std::size_t>();

    // Prepare chunk buffer: need 3-dim array storing 2-dim vectors
    constexpr std::size_t vec_dims = 2;     // we only need to store E_r and E_z
    using chunk_t = typename SpatialSymmetry::Cylindrical<scalar_t, vec_dims>::chunk_t;
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t, vec_dims>::darr_t;
    
    RZTVector<std::size_t> chunk_size(max_pts_in_chunk);
    chunk_t chunk_buffer(chunk_size);
    
    // Prepare distributed array
    std::size_t cache_depth = 1;   // all chunks are prepared as final, no need for a large cache here
    RZTVector<std::size_t> init_cache_el_shape = chunk_size;
    RZTVector<std::size_t> streamer_chunk_shape(stor::INFTY); streamer_chunk_shape[0] = 1; // serialize one outermost slice at a time
    darr_t darr(gf_path, cache_depth, init_cache_el_shape, streamer_chunk_shape);

    // Loop over chunks, fill them, and register them in the distributed array
    auto fill_and_register_chunk = [&](const RZTIndexVector& chunk_start_ind, const RZTIndexVector& chunk_end_ind){

      // Prepare the start- and end coordinates of this chunk ...
      RZTVector<std::size_t> cur_chunk_size = chunk_end_ind - chunk_start_ind;
      RZTCoordVector chunk_start_coords = start_coords_rzt + chunk_start_ind.template as_type<scalar_t>() * stepsize;
      RZTCoordVector chunk_end_coords = start_coords_rzt + chunk_end_ind.template as_type<scalar_t>() * stepsize;

      // ... and make sure the chunk buffer matches the required shape
      chunk_buffer.resize(cur_chunk_size);

      // Fill the buffer ...
      
      std::cout << "filling chunk with start = " << chunk_start_ind << " and end = " << chunk_end_ind << std::endl;
      std::cout << "chunk_start_coords = " << chunk_start_coords << " chunk_end_coords = " << chunk_end_coords << std::endl;

      // ... and register it
      darr.RegisterChunk(chunk_buffer, chunk_start_ind);
    };    
    IteratorUtils::index_loop_over_chunks(start_inds, end_inds, chunk_size, fill_and_register_chunk);

    // Create the actual Green's function from the sampled data
    CylindricalGreensFunction(start_coords_rzt, end_coords_rzt, stepsize, std::move(darr));
  }
}
