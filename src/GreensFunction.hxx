#include <cmath>
#include "GreensFunction.hh"
#include "Serialization.hh"
#include "Interpolation.hh"
#include "MemoryUtils.hh"

CylindricalGreensFunctionMetadata::CylindricalGreensFunctionMetadata() :
  start_pos_rzt(0), end_pos_rzt(0), sample_interval_rzt(0) { }

CylindricalGreensFunctionMetadata::CylindricalGreensFunctionMetadata(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval) :
  start_pos_rzt(start_pos), end_pos_rzt(end_pos), sample_interval_rzt(sample_interval),
  start_pos_rz{start_pos.r(), start_pos.z()}, end_pos_rz{end_pos.r(), end_pos.z()}, sample_interval_rz{sample_interval.r(), sample_interval.z()} { }

namespace stor {

  template <>
  struct Traits<CylindricalGreensFunctionMetadata> {
    using type = CylindricalGreensFunctionMetadata;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<RZTCoordVector>::serialize(stream, val.start_pos_rzt);
      Traits<RZTCoordVector>::serialize(stream, val.end_pos_rzt);
      Traits<RZTCoordVector>::serialize(stream, val.sample_interval_rzt);
    }

    static type deserialize(std::iostream& stream) {
      RZTCoordVector start_pos_rzt = Traits<RZTCoordVector>::deserialize(stream);
      RZTCoordVector end_pos_rzt = Traits<RZTCoordVector>::deserialize(stream);
      RZTCoordVector sample_interval_rzt = Traits<RZTCoordVector>::deserialize(stream);      
      return CylindricalGreensFunctionMetadata(start_pos_rzt, end_pos_rzt, sample_interval_rzt);
    }
  };  
}

// ---------

CylindricalGreensFunction::CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size) :
  lib_t(path, cache_size), m_meta(), m_meta_path(path / m_meta_filename) {
  load_metadata();   // Load metadata from disk
}

CylindricalGreensFunction::CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
						     Distributed_RZT_ErEz_Array&& data) :
  lib_t(move_path_from(std::move(data))), m_meta(start_pos, end_pos, sample_interval), m_meta_path(GetLibdir() / m_meta_filename) {
  save_metadata();   // Dump metadata to disk right away
}

template <class KernelT, class QuadratureT>
void CylindricalGreensFunction::apply_accumulate(const LineCurrentSegment& seg, scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples,
						 std::vector<scalar_t>& signal) {

  // Max. integration step size delta_t_p
  scalar_t max_itgr_step = 1.0f;
  
  // Determine total number of quadrature intervals (always integer) ...
  const std::size_t num_quadrature_intervals = std::ceil((seg.end_time - seg.start_time) / max_itgr_step);

  // ... and the actual integration step size
  scalar_t itgr_step = (seg.end_time - seg.start_time) / num_quadrature_intervals;

  // The set of integration steps includes the first and last point of the full interval
  const std::size_t num_itgr_steps = num_quadrature_intervals + 1;

  // calculate velocity vector and step size along the trajectory
  XYZCoordVector seg_vel = (seg.end_pos - seg.start_pos) / (seg.end_time - seg.start_time);
  XYZCoordVector seg_step = seg_vel * itgr_step;

  // the current represented by this segment
  XYZCoordVector source_xyz = seg_vel * seg.charge;
  
  // Guess a good integration block size: this is purely for reasons of efficiency and will not change the result
  
  // For the t_sig direction, can just take the average chunk size along t -> convert into number of samples by dividing by t_sig_samp
  // For the t_p direction, take min[average_chunk_size_rz / velocity_rz] -> convert into number of integration steps by dividing by itgr_step

  const std::size_t max_sample_block_size = 10; // number of signal samples in block
  const std::size_t max_itgr_block_size = 10; // number of integration steps in block
  
  NDVecArray<scalar_t, 1, vec_dims> coords_rz(max_itgr_block_size);
  NDVecArray<scalar_t, 1, vec_dims> source_rz(max_itgr_block_size);  
  std::vector<scalar_t, no_init_alloc<scalar_t>> quadrature_weights(max_itgr_block_size);

  std::cout << "HHHHH" << std::endl;
  std::cout << "calculating with num_samples = " << num_samples << std::endl;
  std::cout << "using num_itgr_steps = " << num_itgr_steps << std::endl;
  std::cout << "HHHHH" << std::endl;
  
  // Iterate over (integration, output sample) blocks
  for(std::size_t itgr_block_start = 0; itgr_block_start < num_itgr_steps; itgr_block_start += max_itgr_block_size) {
    std::size_t itgr_block_size = std::min(max_itgr_block_size, num_itgr_steps - itgr_block_start);  // actual size of this integration block
    
    // Precompute a few things that we can then reuse for all sample blocks
    for(std::size_t i_pt = 0; i_pt < itgr_block_size; i_pt++) {

      // Current position along the segment
      std::size_t itgr_pt = i_pt + itgr_block_start;
      XYZCoordVector seg_pos = seg.start_pos + itgr_pt * seg_step;
      
      // precalculate (r, z) coordinates of the segment evaluation points in this integration block
      coord_cart_to_cyl(seg_pos, coords_rz[i_pt]);
      
      // precalculate (r, z) field components of the current source in this integration block
      field_cart_to_cyl(source_xyz, seg_pos, source_rz[i_pt]);
    }

    // Precompute quadrature weights for the evaluation points in this integration block
    quadrature_weights.resize(itgr_block_size);
    QuadratureT::fill_weights(itgr_block_start, itgr_block_start + itgr_block_size, num_itgr_steps, quadrature_weights.begin(), quadrature_weights.end());
    
    // Iterate over sample blocks
    for(std::size_t sample_block_start = 0; sample_block_start < num_samples; sample_block_start += max_sample_block_size) {
      std::size_t sample_block_end = std::min(sample_block_start + max_sample_block_size, num_samples);

      std::cout << "---------" << std::endl;
      std::cout << "now on integration block with start = " << itgr_block_start << ", size = " << itgr_block_size << std::endl;
      std::cout << "now on sample block with start = " << sample_block_start << ", end = " << sample_block_end << std::endl;
      
      // iterate over the segment evaluation points in this integration block
      for(std::size_t i_pt = 0; i_pt < itgr_block_size; i_pt++) {

	// Global index of the current integration point
	std::size_t itgr_pt = i_pt + itgr_block_start;

	// Integration time
	scalar_t t_p = itgr_pt * itgr_step + seg.start_time;
	
	// Find index of the first nonzero output sample that we can determine, i.e. the first sample for which t - t' > 0 in the convolution integral
	std::size_t sample_ind_causal = std::max<int>(0, std::ceil((t_p - t_sig_start) / t_sig_samp));
	std::size_t block_sample_ind_start = std::max(sample_ind_causal, sample_block_start);

	// Number of samples we're computing now
	std::size_t block_num_samples = sample_block_end - block_sample_ind_start;	

	// Time of the first output sample that we are computing in this block
	scalar_t block_t_sig_start = block_sample_ind_start * t_sig_samp + t_sig_start;
	
	// start time in integration block
	scalar_t convolution_t_start = block_t_sig_start - t_p;

	std::cout << " . . . . . . . . " << std::endl;
	std::cout << "start to call accumulate_inner_product" << std::endl;
	std::cout << "t_p = " << t_p << std::endl;
	std::cout << "t_sig_start = " << t_sig_start << std::endl;
	std::cout << "block_t_sig_start = " << block_t_sig_start << std::endl;
	std::cout << "convolution_t_start = " << convolution_t_start << std::endl;
	std::cout << "output sample offset = " << block_sample_ind_start << std::endl;
	std::cout << "coords_rz = " << coords_rz[i_pt] << std::endl;
	std::cout << "t_sig_samp = " << t_sig_samp << std::endl;
	std::cout << "num_samples = " << block_num_samples << std::endl;
	std::cout << "source_rz = " << source_rz[i_pt] << std::endl;
	std::cout << "quadrature_weight = " << quadrature_weights[i_pt] << std::endl;

	// The above should make sure we're only interested in the causal part of the Green's function
	assert(convolution_t_start >= 0.0);
	
	auto block_result = signal.begin() + block_sample_ind_start;	
	accumulate_inner_product<KernelT>(coords_rz[i_pt], convolution_t_start, t_sig_samp, block_num_samples, source_rz[i_pt], block_result,
					  quadrature_weights[i_pt] * itgr_step);

	std::cout << " . . . . . . . . " << std::endl;
      }

      std::cout << "-----------" << std::endl;
    }         
  }  
}

template <class KernelT>
void CylindricalGreensFunction::accumulate_inner_product(const RZCoordVectorView rz_coords, scalar_t t_start, scalar_t t_samp, std::size_t num_samples,
							 const RZFieldVectorView source, std::vector<scalar_t>::iterator result, scalar_t weight) {

  // Buffer to hold the interpolated values
  constexpr std::size_t init_interp_buffer_len = 100;
  NDVecArray<scalar_t, 1, vec_dims> interp_buffer(init_interp_buffer_len);

  // Convert everything from coordinates to floating-point array indices (`f_ind`)
  RZVector rz_f_ind = coords_to_index(rz_coords);
  scalar_t t_start_f_ind = (t_start - m_meta.start_pos_rzt.t()) / m_meta.sample_interval_rzt.t();
  scalar_t t_samp_f_ind = t_samp / m_meta.sample_interval_rzt.t();
  
  // Iterate over all chunks that are required to cover the range [t_start_f_ind, t_end_f_ind)
  RZTVector cur_f_ind(rz_f_ind, t_start_f_ind);

  std::size_t cur_sample_ind = 0;
  while(cur_sample_ind < num_samples) {

    // Time coordinate of the current sample
    cur_f_ind.t() = t_start_f_ind + cur_sample_ind * t_samp_f_ind;
    
    // Fetch the chunk that contains this current location ...
    RZTIndexVector cur_ind = cur_f_ind.template as_type<std::size_t>();
    const metadata_t& meta = m_index.GetChunk(cur_ind);

    // TODO; if meta.chunk_type == ChunkType.all_null, simply skip everything else and jump to the next chunk
    
    const chunk_t& chunk = m_cache.RetrieveChunk(meta);

    // Ensure that the overlap on the loaded chunk is large enough for the kernel that we use
    assert(meta.overlap >= KernelT::support);
    
    // Check how many of the next samples are determined by the currently-loaded chunk
    std::size_t chunk_num_samples = (std::size_t)((RZTIndexVector::t(meta.end_ind) - cur_f_ind.t()) / t_samp_f_ind) + 1;

    // std::cout << "chunk_num_samples = " << chunk_num_samples << std::endl;
    
    // Make sure to not request more samples than needed
    std::size_t samples_to_request = std::min(chunk_num_samples, num_samples - cur_sample_ind);

    // Ensure that we're not going too far (interpolation end point `cur_t_end_f_ind` is exclusive,
    // i.e. can reach up to and including the `end_ind` of the chunk)
    assert(cur_f_ind.t() + (samples_to_request - 1) * t_samp_f_ind <= RZTIndexVector::t(meta.end_ind));
    
    // std::cout << "HHHHHHHHH" << std::endl;
    // std::cout << meta << std::endl;
    // std::cout << "requesting interpolation from start = " << cur_f_ind.t() << " to end = " << cur_f_ind.t() + (samples_to_request - 1) * t_samp_f_ind <<
    //   "; " << samples_to_request << " samples" << std::endl;
    // std::cout << "HHHHHHHHH" << std::endl;
       
    // Convert to chunk-local coordinates
    // TODO: this is very clunky, make better
    RZTIndexVector rzt_chunk_ind_offset{RZTIndexVector::r(meta.loc_ind_offset), RZTIndexVector::z(meta.loc_ind_offset), RZTIndexVector::t(meta.loc_ind_offset)};
    RZTVector rzt_chunk_local_ind = (cur_ind - rzt_chunk_ind_offset).template as_type<scalar_t>() + (cur_f_ind - cur_ind.template as_type<scalar_t>());
    RZCoordVector rz_chunk_local_ind{rzt_chunk_local_ind.r(), rzt_chunk_local_ind.z()};

    // Reset interpolation buffer
    interp_buffer.clear();
    
    // Perform interpolation
    Interpolation::interpolate<KernelT>(chunk, interp_buffer,
					rz_chunk_local_ind, rzt_chunk_local_ind.t(),
					t_samp_f_ind, samples_to_request);

    // Compute inner products and accumulate in output range
    for(std::size_t i = 0; i < samples_to_request; i++) {
      *result += inner_product(interp_buffer[i], source) * weight;
      std::advance(result, 1);
    }
    
    cur_sample_ind += samples_to_request;
  }
}

void CylindricalGreensFunction::coord_cart_to_cyl(const XYZCoordVector& coords_cart, RZCoordVectorView coords_cyl) {
  scalar_t x = coords_cart.x();
  scalar_t y = coords_cart.y();
  RZCoordVectorView::r(coords_cyl) = std::sqrt(x * x + y * y);
  RZCoordVectorView::z(coords_cyl) = coords_cart.z();
}

void CylindricalGreensFunction::field_cart_to_cyl(const XYZFieldVector& field_cart, const XYZCoordVector& coords_cart, RZFieldVectorView field_cyl) {
  scalar_t x = coords_cart.x();
  scalar_t y = coords_cart.y();  
  scalar_t r_xy = std::sqrt(x * x + y * y);

  scalar_t cos_phi = 0.0, sin_phi = 0.0;
  if(r_xy > 0) {
    cos_phi = x / r_xy;
    sin_phi = y / r_xy;
  }
    
  RZFieldVectorView::r(field_cyl) = field_cart.x() * cos_phi + field_cart.y() * sin_phi;
  RZFieldVectorView::z(field_cyl) = field_cart.z();
}

scalar_t CylindricalGreensFunction::inner_product(const RZFieldVectorView field, const RZFieldVectorView source) {
  return field[0] * source[0] + field[1] * source[1];
}

RZCoordVector CylindricalGreensFunction::coords_to_index(const RZCoordVectorView rz_coords) {
  
  RZCoordVector coord_vec(rz_coords);
  return (coord_vec - m_meta.start_pos_rz) / m_meta.sample_interval_rz;
}

RZTCoordVector CylindricalGreensFunction::coords_to_index(const RZTCoordVector& rzt_coords) {
  return (rzt_coords - m_meta.start_pos_rzt) / m_meta.sample_interval_rzt;
}

void CylindricalGreensFunction::save_metadata() {
  std::fstream ofs;
  ofs.open(m_meta_path, std::ios::out | std::ios::binary);
  stor::Traits<CylindricalGreensFunctionMetadata>::serialize(ofs, m_meta);
  ofs.close();
}

void CylindricalGreensFunction::load_metadata() {  
  std::fstream ifs;
  ifs.open(m_meta_path, std::ios::in | std::ios::binary);
  m_meta = stor::Traits<CylindricalGreensFunctionMetadata>::deserialize(ifs);
  ifs.close();
}

std::filesystem::path CylindricalGreensFunction::move_path_from(Distributed_RZT_ErEz_Array&& data) {
  data.Flush(); // ensure that all data is ready to be read from disk
  return data.GetWorkdir();
}
