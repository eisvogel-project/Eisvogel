#include "GreensFunction.hh"
#include "Serialization.hh"
#include "Interpolation.hh"

CylindricalGreensFunctionMetadata::CylindricalGreensFunctionMetadata() :
  start_pos(0), end_pos(0), sample_interval(0) { }

CylindricalGreensFunctionMetadata::CylindricalGreensFunctionMetadata(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval) :
  start_pos(start_pos), end_pos(end_pos), sample_interval(sample_interval) { }

namespace stor {

  template <>
  struct Traits<CylindricalGreensFunctionMetadata> {
    using type = CylindricalGreensFunctionMetadata;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<RZTCoordVector>::serialize(stream, val.start_pos);
      Traits<RZTCoordVector>::serialize(stream, val.end_pos);
      Traits<RZTCoordVector>::serialize(stream, val.sample_interval);
    }

    static type deserialize(std::iostream& stream) {
      RZTCoordVector start_pos = Traits<RZTCoordVector>::deserialize(stream);
      RZTCoordVector end_pos = Traits<RZTCoordVector>::deserialize(stream);
      RZTCoordVector sample_interval = Traits<RZTCoordVector>::deserialize(stream);      
      return CylindricalGreensFunctionMetadata(start_pos, end_pos, sample_interval);
    }
  };  
}

// ---------

CylindricalGreensFunction::CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size) :
  m_meta_path(path / m_meta_filename), m_lib(path, cache_size) {
  load_metadata();   // Load metadata from disk
}

CylindricalGreensFunction::CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
						     Distributed_RZT_ErEz_Array&& data) :
  m_meta(start_pos, end_pos, sample_interval), m_lib(move_path_from(std::move(data))) {
  save_metadata();   // Dump metadata to disk right away
}

template <class KernelT>
void CylindricalGreensFunction::accumulate_inner_product(const XYZCoordVector& coords, scalar_t t_start, scalar_t t_end, scalar_t t_samp,
							 const XYZTFieldVector& current, std::vector<scalar_t>::iterator result) {

  // Convert (x, y, z) coordinates into indices
  
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
