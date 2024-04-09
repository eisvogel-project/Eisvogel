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
void CylindricalGreensFunction::accumulate_inner_product(const RZCoordVector& rz_coords, scalar_t t_start, scalar_t t_end, scalar_t t_samp,
							 const XYZFieldVector& current, std::vector<scalar_t>::iterator result) {

  // Convert coordinates to index values
  RZTVector rzt_start_coords(rz_coords, t_start);  
  RZTVector rzt_start_inds = coords_to_index(rzt_start_coords);
  
  RZVector rz_inds = coords_to_index(rz_coords);
    
  
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
