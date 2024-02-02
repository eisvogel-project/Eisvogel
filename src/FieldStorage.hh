#ifndef __FIELD_STORAGE_HH
#define __FIELD_STORAGE_HH

#include <string>
#include <filesystem>
#include "Eisvogel/DistributedNDArray.hh"

class FieldStorage {
  
public:
  virtual IndexVector shape() = 0;
};

class RZFieldStorage : public FieldStorage {

public:
  using storage_t = DistributedNDArray<scalar_t, 3>;

public:

  RZFieldStorage(std::filesystem::path storage_path, std::size_t cache_size);
  
  scalar_t E_r(IndexVector& ind);
  scalar_t E_z(IndexVector& ind);  
  scalar_t E_phi(IndexVector& ind);

  IndexVector shape();

  // vector E_rzphi(IndexVector& ind);
  
  // calculate from the stored components
  // scalar E_x(IndexVector& ind);
  // scalar E_y(IndexVector& ind);
  // scalar E_z(IndexVector& ind);

  // vector E_xyz(IndexVector& ind);
  
private:
  
  std::shared_ptr<storage_t> m_E_r;
  std::shared_ptr<storage_t> m_E_z;  
};

#endif
