#ifndef __FIELD_STORAGE_HH
#define __FIELD_STORAGE_HH

#include <string>
#include <filesystem>
#include "DistributedNDArray.hh"

class FieldStorage {
  
public:
  virtual IndexVector shape() = 0;
};

class RZFieldStorage : public FieldStorage {

public:
  using serializer_t = stor::DefaultSerializer;
  using storage_t = DistributedScalarNDArray<scalar_t, 3, serializer_t>;
  
public:

  RZFieldStorage(std::filesystem::path storage_path, std::size_t cache_size);
  
  scalar_t E_r(IndexVector& ind);
  scalar_t E_z(IndexVector& ind);  
  scalar_t E_phi(IndexVector& ind);

  IndexVector shape();

  void MakeIndexPersistent();

  template <class ChunkT>
  void RegisterChunk(const ChunkT& chunk_E_r, const ChunkT& chunk_E_z, const IndexVector chunk_start_inds) {
    m_E_r -> RegisterChunk(chunk_E_r, chunk_start_inds);
    m_E_z -> RegisterChunk(chunk_E_z, chunk_start_inds);
  };

  void RebuildChunks(const IndexVector& requested_chunk_size);
  void MergeChunks(std::size_t dim_to_merge, std::size_t max_dimsize);
  
  // vector E_rzphi(IndexVector& ind);
  
  // calculate from the stored components
  // scalar E_x(IndexVector& ind);
  // scalar E_y(IndexVector& ind);
  // scalar E_z(IndexVector& ind);

  // vector E_xyz(IndexVector& ind);
  
private:
  
  std::shared_ptr<storage_t> m_E_r;
  std::shared_ptr<storage_t> m_E_z;

  std::shared_ptr<serializer_t> m_ser;
};

#endif
