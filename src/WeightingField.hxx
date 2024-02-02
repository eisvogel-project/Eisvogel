#include "Eisvogel/Serialization.hh"

template <typename SymmetryT>
TRZFieldIndexer<SymmetryT>::TRZFieldIndexer(std::filesystem::path index_path) {

  std::fstream ifs;
  std::filesystem::path meta_path = index_path / "wf_meta.bin";
  ifs.open(meta_path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);
  m_start_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());
  m_end_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());

}

template <typename SymmetryT>
void TRZFieldIndexer<SymmetryT>::SetShape(IndexVector shape) {
  m_shape = std::make_shared<IndexVector>(shape);
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
  
  IndexVector ind = {(std::size_t)ind_t, (std::size_t)ind_z, (std::size_t)ind_r};	
  return ind;  
}

// ---------------

template <class FieldIndexerT, class FieldStorageT>
WeightingField<FieldIndexerT, FieldStorageT>::WeightingField(std::string wf_path) : m_wf_path(wf_path) {

  m_field_indexer = std::make_shared<FieldIndexerT>(wf_path);  
  m_field_storage = std::make_shared<FieldStorageT>(wf_path, 10);

  m_field_indexer->SetShape(m_field_storage->shape());
}
