#include "GreensFunction.hh"
#include "Serialization.hh"

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

template <class KernelT>
CylindricalGreensFunction<KernelT>::CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size) :
  m_meta_path(path / m_meta_filename), m_lib(path, cache_size) {

  // Load metadata from disk
  load_metadata();
}

template <class KernelT>
CylindricalGreensFunction<KernelT>::CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
							      Distributed_RZT_ErEz_Array&& data) :
  m_meta(start_pos, end_pos, sample_interval), m_lib(move_path_from(std::move(data))) {

  // Dump metadata to disk right away
  save_metadata();
}

template <class KernelT>
void CylindricalGreensFunction<KernelT>::accumulate_inner(const XYZCoordVector& coords, scalar_t t_start, scalar_t t_end, scalar_t t_samp,
							  const XYZTFieldVector& current, std::vector<scalar_t>::iterator result) {

}

template <class KernelT>
void CylindricalGreensFunction<KernelT>::save_metadata() {
  std::fstream ofs;
  ofs.open(m_meta_path, std::ios::out | std::ios::binary);
  stor::Traits<CylindricalGreensFunctionMetadata>::serialize(ofs, m_meta);
  ofs.close();
}

template <class KernelT>
void CylindricalGreensFunction<KernelT>::load_metadata() {  
  std::fstream ifs;
  ifs.open(m_meta_path, std::ios::in | std::ios::binary);
  m_meta = stor::Traits<CylindricalGreensFunctionMetadata>::deserialize(ifs);
  ifs.close();
}

template <class KernelT>
std::filesystem::path CylindricalGreensFunction<KernelT>::move_path_from(Distributed_RZT_ErEz_Array&& data) {
  data.Flush(); // ensure that all data is ready to be read from disk
  return data.GetWorkdir();
}
