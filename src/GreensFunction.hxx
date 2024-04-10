#include "GreensFunction.hh"
#include "Serialization.hh"
#include "Interpolation.hh"

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

template <class KernelT>
void CylindricalGreensFunction::accumulate_inner_product(const RZCoordVector& rz_coords, scalar_t t_start, scalar_t t_samp, std::size_t num_samples,
							 const RZFieldVector& source, std::vector<scalar_t>::iterator result) {

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

    // Reset interpolation buffer
    interp_buffer.clear();
    
    // Perform interpolation
    // TODO: the index conversion is very clunky ... let's find a nice way to do this
    RZCoordVector rz_chunk_ind_offset{(scalar_t)RZTIndexVector::r(meta.loc_ind_offset), (scalar_t)RZTIndexVector::z(meta.loc_ind_offset)};
    Interpolation::interpolate<KernelT>(chunk, interp_buffer,
					rz_f_ind - rz_chunk_ind_offset,  // convert to chunk-local coordinates
					cur_f_ind.t() - RZTIndexVector::t(meta.loc_ind_offset),  // convert to chunk-local coordinates
					t_samp_f_ind, samples_to_request);

    // Compute inner products and accumulate in output range
    for(std::size_t i = 0; i < samples_to_request; i++) {
      *result += inner_product(interp_buffer[i], source);
      std::advance(result, 1);
    }
    
    cur_sample_ind += samples_to_request;
  }
}

scalar_t CylindricalGreensFunction::inner_product(const view_t& field, const RZFieldVector& source) {
  return field[0] * source[0] + field[1] * source[1];
}

RZCoordVector CylindricalGreensFunction::coords_to_index(const RZCoordVector& coords) {
  return (coords - m_meta.start_pos_rz) / m_meta.sample_interval_rz;
}

RZTCoordVector CylindricalGreensFunction::coords_to_index(const RZTCoordVector& coords) {
  return (coords - m_meta.start_pos_rzt) / m_meta.sample_interval_rzt;
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
