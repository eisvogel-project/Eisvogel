#include <cassert>
#include "IteratorUtils.hh"
#include "DistributedNDVecArray.hh"
#include "GreensFunction.hh"
#include "MathUtils.hh"
#include "Symmetry.hh"

namespace {

  void EvaluateElectricDipoleGreensFunction(
      const RZTCoordVector& start_coords, const RZTCoordVector& stepsize,
      scalar_t ior, scalar_t filter_t_peak, unsigned int filter_order,
      scalar_t r_min, SpatialSymmetry::Cylindrical<scalar_t>::chunk_t& buffer)
  {

    scalar_t Qw = 1.0;
    scalar_t eps0 = 1.0; // vacuum dielectric constant
    scalar_t ds = 1.0;
    scalar_t c = 1.0; // speed of light in vacuum

    auto filtered_theta = [&](scalar_t t) -> scalar_t {
      if (t <= 0) {
        return 0.0;
      }
      return 1.0 - MathUtils::incomplete_gamma(
                       1 + filter_order, filter_order * t / filter_t_peak) /
                       std::exp(std::lgamma(filter_order + 1));
    };

    auto filtered_delta = [&](scalar_t t) -> scalar_t {
      if (t <= 0) {
        return 0.0;
      }
      return std::pow(t / filter_t_peak * filter_order, filter_order) *
             std::exp(-t / filter_t_peak * filter_order) /
             (filter_t_peak * std::exp(std::lgamma(filter_order)));
    };

    auto filtered_delta_prime = [&](scalar_t t) -> scalar_t {
      if (t <= 0) {
        return 0.0;
      }
      return filtered_delta(t) * (filter_t_peak - t) * filter_order /
             (filter_t_peak * t);
    };

    // Weighting field in spherical coordinates
    auto E_r = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(r_xy * r_xy + z * z);
      if (r < r_min) {
        return std::nan("");
      }
      scalar_t t_prop = r * ior / c, t_del = t - t_prop;
      scalar_t cos_theta = z / r;

      return -2.0 * Qw * ds / (eps0 * 4 * M_PI) * cos_theta / std::pow(r, 3) *
             (filtered_theta(t_del) + t_prop * filtered_delta(t_del));
    };

    auto E_theta = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(r_xy * r_xy + z * z);
      if (r < r_min) {
        return std::nan("");
      }
      scalar_t t_prop = r * ior / c, t_del = t - t_prop;
      scalar_t sin_theta = r_xy / r;
      return -Qw * ds / (eps0 * 4 * M_PI) * sin_theta / std::pow(r, 3) *
             (filtered_theta(t_del) + t_prop * filtered_delta(t_del) +
              std::pow(t_prop, 2) * filtered_delta_prime(t_del));
    };

    // Weighting field in cylindrical coordinates
    auto E_rxy = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(r_xy * r_xy + z * z);
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * sin_theta + E_theta(t, r_xy, z) * cos_theta;
    };

    auto E_z = [&](scalar_t t, scalar_t r_xy, scalar_t z) -> scalar_t {
      r_xy = std::fabs(r_xy);
      scalar_t r = std::sqrt(r_xy * r_xy + z * z);
      scalar_t cos_theta = z / r, sin_theta = r_xy / r;
      return E_r(t, r_xy, z) * cos_theta - E_theta(t, r_xy, z) * sin_theta;
    };

    using view_t = typename SpatialSymmetry::Cylindrical<scalar_t>::view_t;
    auto evaluator = [&](const RZTIndexVector& ind, view_t cur_element) {
      // Coordinates of current point
      RZTCoordVector cur_pt =
          start_coords + stepsize * (ind.template as_type<scalar_t>());

      // Fields at this location
      scalar_t cur_E_rxy = E_rxy(cur_pt.t(), cur_pt.r(), cur_pt.z());
      scalar_t cur_E_z = E_z(cur_pt.t(), cur_pt.r(), cur_pt.z());

      // Store
      cur_element[0] = cur_E_rxy;
      cur_element[1] = cur_E_z;
    };

    // evaluate the Green's function and fill the buffer
    buffer.index_loop_over_elements(evaluator);
  }
} // namespace

namespace GreensFunctionCalculator::Analytic {

  void ElectricDipole(std::filesystem::path gf_path,
                      const RZCoordVector& start_coords,
                      const RZCoordVector& end_coords, scalar_t t_end,
                      scalar_t ior, scalar_t filter_t_peak,
                      unsigned int filter_order, scalar_t r_min,
                      scalar_t os_factor, std::size_t max_pts_in_chunk)
  {

    // make sure to start from scratch
    std::filesystem::remove_all(gf_path);

    // compute required step size for sampling of weighting field
    scalar_t c = 1.0; // speed of light in vacuum
    scalar_t fmax = (scalar_t)(filter_order) / (2 * M_PI * filter_t_peak) *
                    std::sqrt(std::pow(2.0, 1.0 / (filter_order + 1)) - 1);
    scalar_t lambda_min = c / (fmax * ior);
    scalar_t delta_t = 1.0 / (2 * fmax * os_factor);
    scalar_t delta_pos = lambda_min / (2.0 * os_factor);

    scalar_t t_start = 0.0;

    std::cout << "---------------------------" << std::endl;
    std::cout << "Using oversampling factor = " << os_factor << std::endl;
    std::cout << "delta_t = " << delta_t << std::endl;
    std::cout << "delta_r = delta_z = " << delta_pos << std::endl;
    std::cout << "start_coords: t = " << t_start << ", r = " << start_coords.r()
              << ", z = " << start_coords.z() << std::endl;
    std::cout << "end_coords: t = " << t_end << ", r = " << end_coords.r()
              << ", z = " << end_coords.z() << std::endl;
    std::cout << "---------------------------" << std::endl;

    RZTCoordVector start_coords_rzt{start_coords.r(), start_coords.z(),
                                    t_start};
    RZTCoordVector end_coords_rzt{end_coords.r(), end_coords.z(), t_end};

    // Determine step size so that an integer number of samples fits into the
    // domain of the Green's function
    RZTVector<scalar_t> stepsize_requested{delta_pos, delta_pos, delta_t};
    RZTVector<std::size_t> number_pts =
        ((end_coords_rzt - start_coords_rzt) / stepsize_requested)
            .template as_type<std::size_t>();
    RZTVector<scalar_t> stepsize = (end_coords_rzt - start_coords_rzt) /
                                   number_pts.template as_type<scalar_t>();

    // Resulting start- and end indices
    RZTIndexVector start_inds(0);
    RZTIndexVector end_inds = number_pts + 1;

    // Prepare chunk buffer: need 3-dim array storing 2-dim vectors
    using chunk_t = typename SpatialSymmetry::Cylindrical<scalar_t>::chunk_t;
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;

    RZTVector<std::size_t> chunk_size(max_pts_in_chunk);
    chunk_t chunk_buffer(chunk_size);

    // Prepare distributed array
    std::size_t cache_depth =
        5; // all chunks are prepared as final, no need for a large cache here
    RZTVector<std::size_t> init_cache_el_shape = chunk_size;
    RZTVector<std::size_t> streamer_chunk_shape(stor::INFTY);
    streamer_chunk_shape[0] = 1; // serialize one outermost slice at a time
    darr_t darr(gf_path, cache_depth, init_cache_el_shape,
                streamer_chunk_shape);

    // Loop over chunks, fill them, and register them in the distributed array
    auto fill_and_register_chunk = [&](const RZTIndexVector& chunk_start_ind,
                                       const RZTIndexVector& chunk_end_ind) {
      // Prepare the start- and end coordinates of this chunk ...
      RZTVector<std::size_t> cur_chunk_size = chunk_end_ind - chunk_start_ind;
      RZTCoordVector chunk_start_coords =
          start_coords_rzt +
          chunk_start_ind.template as_type<scalar_t>() * stepsize;

      // ... and make sure the chunk buffer matches the required shape
      chunk_buffer.resize(cur_chunk_size);

      // Fill the buffer ...
      EvaluateElectricDipoleGreensFunction(chunk_start_coords, stepsize, ior,
                                           filter_t_peak, filter_order, r_min,
                                           chunk_buffer);

      // ... and register it
      darr.RegisterChunk(chunk_buffer, chunk_start_ind);
    };
    IteratorUtils::index_loop_over_chunks(start_inds, end_inds, chunk_size,
                                          fill_and_register_chunk);

    // Reshape the chunks to make sure the overlap required for the
    // interpolation is respected
    std::size_t overlap = 2;
    std::filesystem::path workdir_tmp = "./darr_test_tmp";
    darr.RebuildChunks(
        chunk_size, workdir_tmp, overlap,
        SpatialSymmetry::Cylindrical<scalar_t>::boundary_evaluator);

    // Create the actual Green's function from the sampled data
    CylindricalGreensFunction(start_coords_rzt, end_coords_rzt, stepsize,
                              std::move(darr));
  }
} // namespace GreensFunctionCalculator::Analytic
