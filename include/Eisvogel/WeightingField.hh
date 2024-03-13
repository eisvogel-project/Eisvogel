#ifndef __WEIGHTING_FIELD_HH
#define __WEIGHTING_FIELD_HH

#include <memory>
#include <string>
#include <filesystem>

#include "Common.hh"
#include "CoordUtils.hh"
#include "Interpolator.hh"
#include "Kernels.hh"

#include "FieldStorage.hh"

template <typename SymmetryT>
class TRZFieldIndexer {

public:

  TRZFieldIndexer(std::filesystem::path index_path, FieldStorage& storage);
  TRZFieldIndexer(std::filesystem::path index_path, const CoordVector& start_coords, const CoordVector& end_coords);
  
  CoordVector GetFieldIndexFromCoord(CoordVector pos_txyz);
  IndexVector GetFieldIndex(GridVector inds_trz);

  DeltaVector GetSamplingIntervals() const;
  CoordVector GetStartCoords() const;
  CoordVector GetEndCoords() const;

  void MakePersistent();

private:
  
  std::filesystem::path m_meta_path;
  
  std::shared_ptr<IndexVector> m_shape;
  std::shared_ptr<CoordVector> m_start_coords;
  std::shared_ptr<CoordVector> m_end_coords;  
};

template <class FieldIndexerT, class FieldStorageT>
class WeightingField {

public:

  WeightingField(std::string wf_path);
  WeightingField(std::string wf_path, const CoordVector& start_coords, const CoordVector& end_coords);

  template <typename ...Params>
  void RegisterChunk(Params&&... params);

  void MakeMetadataPersistent();
  
  // Accessors for field components
  template <typename KernelT = DefaultKernel>
  scalar_t E_r(CoordVector pos);

  template <typename KernelT = DefaultKernel>
  scalar_t E_z(CoordVector pos);

  template <typename KernelT = DefaultKernel>
  scalar_t E_phi(CoordVector pos);

  // // Vectorized accessors for field components
  // template <typename KernelT = DefaultKernel>
  // std::vector<scalar_t> E_r(std::vector<CoordVector>& pos);

  // template <typename KernelT = DefaultKernel>
  // std::vector<scalar_t> E_z(std::vector<CoordVector>& pos);

  // template <typename KernelT = DefaultKernel>
  // std::vector<scalar_t> E_phi(std::vector<CoordVector>& pos);
  
  // template <typename KernelT = DefaultKernel>
  // scalar_t Eip(CoordVector& pos, CoordVector& vec);

  // template <typename KernelT = DefaultKernel>
  // std::vector<scalar_t> Eip(std::vector<CoordVector>& pos, std::vector<CoordVector>& vec);
  
  DeltaVector GetSamplingIntervals() const;
  CoordVector GetStartCoords() const;
  scalar_t GetStartCoords(std::size_t dim) const;
  CoordVector GetEndCoords() const;
  scalar_t GetEndCoords(std::size_t dim) const;

  void RebuildChunks(const IndexVector& requested_chunk_size);
  void MergeChunks(std::size_t dim_to_merge, std::size_t max_dimsize);
  
private:
  
  template <typename KernelT, typename GetterT,
	    typename ValueT = std::invoke_result_t<GetterT, IndexVector&>>
  ValueT eval(CoordVector pos, GetterT getter);
  
private:

  std::string m_wf_path;
  
  std::shared_ptr<FieldStorageT> m_field_storage;
  std::shared_ptr<FieldIndexerT> m_field_indexer;  
};

#include "WeightingField.hxx"

#include "Symmetry.hh"
using CylindricalWeightingField = WeightingField<TRZFieldIndexer<CylindricalSymmetry>, RZFieldStorage>;

#endif
