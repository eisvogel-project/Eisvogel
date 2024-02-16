#include "FieldStorage.hh"

RZFieldStorage::RZFieldStorage(std::filesystem::path storage_path, std::size_t cache_size) {

  std::filesystem::path storage_path_E_r = storage_path / "E_r";
  std::filesystem::path storage_path_E_z = storage_path / "E_z";

  m_ser = std::make_shared<serializer_t>();
  
  m_E_r = std::make_shared<storage_t>(storage_path_E_r, cache_size, *m_ser);
  m_E_z = std::make_shared<storage_t>(storage_path_E_z, cache_size, *m_ser);
}

IndexVector RZFieldStorage::shape() {
  return m_E_r -> shape();
}

scalar_t RZFieldStorage::E_r(IndexVector& ind) {
  return m_E_r -> operator()(ind);
}

scalar_t RZFieldStorage::E_z(IndexVector& ind) {
  return m_E_z -> operator()(ind);
}

scalar_t RZFieldStorage::E_phi(IndexVector& ind) {
  return 0.0;
}

void RZFieldStorage::MakeIndexPersistent() {
  m_E_r -> MakeIndexPersistent();
  m_E_z -> MakeIndexPersistent();
}

void RZFieldStorage::RebuildChunks(const IndexVector& requested_chunk_size) {
  m_E_r -> RebuildChunks(requested_chunk_size);
  m_E_z -> RebuildChunks(requested_chunk_size);
}

void RZFieldStorage::MergeChunks(std::size_t dim_to_merge, std::size_t max_dimsize) {
  m_E_r -> MergeChunks(dim_to_merge, max_dimsize);
  m_E_z -> MergeChunks(dim_to_merge, max_dimsize);
}

// calculate from the stored components
// scalar E_x(IndexVector& ind) { ... };
// scalar E_y(IndexVector& ind) { ... };
// scalar E_z(IndexVector& ind) { ... };
