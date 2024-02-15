#include "Serialization.hh"

template <typename SymmetryT>
TRZFieldIndexer<SymmetryT>::TRZFieldIndexer(std::filesystem::path index_path, FieldStorage& storage) :
  m_meta_path(index_path / "wf_meta.bin") {

  std::fstream ifs;
  ifs.open(m_meta_path, std::ios::in | std::ios::binary);
  stor::DefaultSerializer iser;
  m_start_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>(ifs));
  m_end_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>(ifs));
  m_shape = std::make_shared<IndexVector>(storage.shape());
}

template <typename SymmetryT>
TRZFieldIndexer<SymmetryT>::TRZFieldIndexer(std::filesystem::path index_path, const CoordVector& start_coords, const CoordVector& end_coords) :
  m_meta_path(index_path / "wf_meta.bin") {

  m_start_coords = std::make_shared<CoordVector>(start_coords);
  m_end_coords = std::make_shared<CoordVector>(end_coords);
}

template <typename SymmetryT>
void TRZFieldIndexer<SymmetryT>::MakePersistent() {
  
  std::fstream ofs;
  ofs.open(m_meta_path, std::ios::out | std::ios::binary);  
  stor::DefaultSerializer oser;
  oser.serialize(ofs, *m_start_coords);
  oser.serialize(ofs, *m_end_coords);
}

template <typename SymmetryT>
CoordVector TRZFieldIndexer<SymmetryT>::GetSamplingIntervals() const {
  return (*m_end_coords - *m_start_coords) / (CoordVector)*m_shape;
}

template <typename SymmetryT>
CoordVector TRZFieldIndexer<SymmetryT>::GetStartCoords() const {
  return *m_start_coords;
}

template <typename SymmetryT>
CoordVector TRZFieldIndexer<SymmetryT>::GetEndCoords() const {
  return *m_end_coords;
}

template <typename SymmetryT>
CoordVector TRZFieldIndexer<SymmetryT>::GetFieldIndexFromCoord(CoordVector pos_txyz) {
  return SymmetryT::GetOrbitIndex(pos_txyz, *m_start_coords, *m_end_coords, *m_shape);
}

template <typename SymmetryT>
IndexVector TRZFieldIndexer<SymmetryT>::GetFieldIndex(GridVector inds_trz) {
  
  int ind_t = inds_trz(0);
  int ind_z = inds_trz(1);
  int ind_r = inds_trz(2);	
  
  if(ind_r < 0) {
    ind_r = -ind_r;
  }
  if(ind_t < 0) {
    ind_t = 0;
  }
  
  IndexVector ind = {(std::size_t)ind_t, (std::size_t)ind_z, (std::size_t)ind_r};	
  return ind;  
}

// ---------------

template <class FieldIndexerT, class FieldStorageT>
WeightingField<FieldIndexerT, FieldStorageT>::WeightingField(std::string wf_path, const CoordVector& start_coords, const CoordVector& end_coords) :
  m_wf_path(wf_path) {

  // we might have to create the directory structure if this is an empty weighting field
  if(!std::filesystem::exists(m_wf_path)) {
    std::filesystem::create_directory(m_wf_path);
  }
  
  m_field_storage = std::make_shared<FieldStorageT>(wf_path, 10);
  m_field_indexer = std::make_shared<FieldIndexerT>(wf_path, start_coords, end_coords);  
}

template <class FieldIndexerT, class FieldStorageT>
WeightingField<FieldIndexerT, FieldStorageT>::WeightingField(std::string wf_path) : m_wf_path(wf_path) {

  // we might have to create the directory structure if this is an empty weighting field
  if(!std::filesystem::exists(m_wf_path)) {
    throw std::runtime_error("No usable weighting field found at this location!");
  }
  
  m_field_storage = std::make_shared<FieldStorageT>(wf_path, 10);
  m_field_indexer = std::make_shared<FieldIndexerT>(wf_path, *m_field_storage);  
}

template <class FieldIndexerT, class FieldStorageT>
template <typename ...Params>
void WeightingField<FieldIndexerT, FieldStorageT>::RegisterChunk(Params&&... params) {
  m_field_storage -> RegisterChunk(std::forward<Params>(params)...);
}

template <class FieldIndexerT, class FieldStorageT>
void WeightingField<FieldIndexerT, FieldStorageT>::MakeMetadataPersistent() {
  m_field_storage -> MakeIndexPersistent();
  m_field_indexer -> MakePersistent();
}

template <class FieldIndexerT, class FieldStorageT>
DeltaVector WeightingField<FieldIndexerT, FieldStorageT>::GetSamplingIntervals() const {
  return m_field_indexer -> GetSamplingIntervals();
}

template <class FieldIndexerT, class FieldStorageT>
CoordVector WeightingField<FieldIndexerT, FieldStorageT>::GetStartCoords() const {
  return m_field_indexer -> GetStartCoords();
}

template <class FieldIndexerT, class FieldStorageT>
scalar_t WeightingField<FieldIndexerT, FieldStorageT>::GetStartCoords(std::size_t dim) const {
  return GetStartCoords()(dim);
}

template <class FieldIndexerT, class FieldStorageT>
CoordVector WeightingField<FieldIndexerT, FieldStorageT>::GetEndCoords() const {
  return m_field_indexer -> GetEndCoords();
}

template <class FieldIndexerT, class FieldStorageT>
scalar_t WeightingField<FieldIndexerT, FieldStorageT>::GetEndCoords(std::size_t dim) const {
  return GetEndCoords()(dim);
}

template <class FieldIndexerT, class FieldStorageT>
void WeightingField<FieldIndexerT, FieldStorageT>::RebuildChunks(const IndexVector& requested_chunk_size) {
  m_field_storage -> RebuildChunks(requested_chunk_size);
}

template <class FieldIndexerT, class FieldStorageT>
template <typename KernelT>
scalar_t WeightingField<FieldIndexerT, FieldStorageT>::E_r(CoordVector pos) {    
  return eval<KernelT>(pos, [&](auto& arg){return m_field_storage -> E_r(arg);});
};

template <class FieldIndexerT, class FieldStorageT>
template <typename KernelT>
scalar_t WeightingField<FieldIndexerT, FieldStorageT>::E_z(CoordVector pos) {    
  return eval<KernelT>(pos, [&](auto& arg){return m_field_storage -> E_z(arg);});
};

template <class FieldIndexerT, class FieldStorageT>
template <typename KernelT>
scalar_t WeightingField<FieldIndexerT, FieldStorageT>::E_phi(CoordVector pos) {    
  return eval<KernelT>(pos, [&](auto& arg){return m_field_storage -> E_phi(arg);});
};

template <class FieldIndexerT, class FieldStorageT>
template <typename KernelT, typename GetterT, typename ValueT>
ValueT WeightingField<FieldIndexerT, FieldStorageT>::eval(CoordVector pos, GetterT getter) {
  
  // Weighting fields are nonzero only for t > 0
  if(CoordUtils::getT(pos) <= 0) {
    return 0.0;
  }
  
  CoordVector frac_inds = m_field_indexer -> GetFieldIndexFromCoord(pos);
  
  auto to_interpolate = [&](GridVector& vec) -> ValueT {
    IndexVector inds = m_field_indexer -> GetFieldIndex(vec);
    return getter(inds);
  };
  
  return InterpolateFunc<KernelT>(to_interpolate, frac_inds);        
}
