#include "Eisvogel/Serialization.hh"

template <typename SymmetryT>
TRZFieldIndexer<SymmetryT>::TRZFieldIndexer(std::filesystem::path index_path, FieldStorage& storage) {

  std::fstream ifs;
  std::filesystem::path meta_path = index_path / "wf_meta.bin";
  ifs.open(meta_path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);
  m_start_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());
  m_end_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());
  m_shape = std::make_shared<IndexVector>(storage.shape());

  m_shape->print();
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

  m_field_storage = std::make_shared<FieldStorageT>(wf_path, 10);
  m_field_indexer = std::make_shared<FieldIndexerT>(wf_path, *m_field_storage);  
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
  
  return InterpolateFuncNew<KernelT>(to_interpolate, frac_inds);        
}
