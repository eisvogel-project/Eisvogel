#ifndef __WEIGHTING_FIELD_HH
#define __WEIGHTING_FIELD_HH

#include <type_traits>
#include "NDArray.hh"
#include "Common.hh"
#include "CoordUtils.hh"

#include <iostream>

template <class StorageT = DenseNDArray<scalar_t, 3>>
class WeightingField {

public:

  using index_t = std::size_t;
  using shape_t = typename StorageT::shape_t;
 
  WeightingField(StorageT&& E_r, StorageT&& E_z, StorageT&& E_phi,
		 CoordVector start_coords, CoordVector end_coords) :
    m_E_r(E_r), m_E_z(E_z), m_E_phi(E_phi), m_start_coords(start_coords), m_end_coords(end_coords) {

    if((m_E_r.shape() == m_E_z.shape()) && (m_E_z.shape() == m_E_phi.shape())) {
      m_shape = m_E_r.shape();
    }
    else {
      throw;
    }
  };

  scalar_t E_r(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_r(ind_t, ind_r, ind_z);
  };

  const StorageT& E_r() const {return m_E_r;};
  const StorageT& E_z() const {return m_E_z;};
  const StorageT& E_phi() const {return m_E_phi;};

  scalar_t E_z(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_z(ind_t, ind_r, ind_z);
  };

  scalar_t E_phi(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_phi(ind_t, ind_r, ind_z);
  };

private:

  shape_t m_shape;

  StorageT m_E_r;
  StorageT m_E_z;
  StorageT m_E_phi;

  CoordVector m_start_coords;
  CoordVector m_end_coords;
};

#endif
