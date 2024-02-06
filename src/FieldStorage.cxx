#include "FieldStorage.hh"

RZFieldStorage::RZFieldStorage(std::filesystem::path storage_path, std::size_t cache_size) {

  std::filesystem::path storage_path_E_r = storage_path / "E_r";
  std::filesystem::path storage_path_E_z = storage_path / "E_z";

  std::shared_ptr<stor::SparseSerializer> ser = std::make_shared<stor::SparseSerializer>();
  
  m_E_r = std::make_shared<storage_t>(storage_path_E_r, cache_size, *ser);
  m_E_z = std::make_shared<storage_t>(storage_path_E_z, cache_size, *ser);
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

void RZFieldStorage::RegisterChunk(const chunk_t& chunk_E_r, const chunk_t& chunk_E_z, const IndexVector chunk_start_inds) {
  m_E_r -> RegisterChunk(chunk_E_r, chunk_start_inds);
  m_E_z -> RegisterChunk(chunk_E_z, chunk_start_inds);
}

void RZFieldStorage::MakeIndexPersistent() {
  m_E_r -> MakeIndexPersistent();
  m_E_z -> MakeIndexPersistent();
}

// calculate from the stored components
// scalar E_x(IndexVector& ind) { ... };
// scalar E_y(IndexVector& ind) { ... };
// scalar E_z(IndexVector& ind) { ... };
