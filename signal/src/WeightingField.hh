#ifndef __WEIGHTING_FIELD_HH
#define __WEIGHTING_FIELD_HH

#include <type_traits>
#include "NDArray.hh"

#include <iostream>

template <template<typename, std::size_t> class StorageT, typename ValueT = float>
class WeightingField {

public:

  using index_t = std::size_t;
  using shape_t = typename StorageT<ValueT, 3>::shape_t;
 
  WeightingField(StorageT<ValueT, 3>&& E_r, StorageT<ValueT, 3>&& E_z, StorageT<ValueT, 3>&& E_phi,
		 ValueT t_min, ValueT t_max, ValueT r_min, ValueT r_max, ValueT z_min, ValueT z_max) :
    m_E_r(E_r), m_E_z(E_z), m_E_phi(E_phi), 
    m_t_min(t_min), m_t_max(t_max), m_r_min(r_min), m_r_max(r_max), m_z_min(z_min), m_z_max(z_max) { 

    if((m_E_r.shape() == m_E_z.shape()) && (m_E_z.shape() == m_E_phi.shape())) {
      m_shape = m_E_r.shape();
    }
    else {
      throw;
    }
  };

  ValueT E_r(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_r(ind_t, ind_r, ind_z);
  };

  ValueT E_z(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_z(ind_t, ind_r, ind_z);
  };

  ValueT E_phi(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_phi(ind_t, ind_r, ind_z);
  };

private:

  shape_t m_shape;

  StorageT<ValueT, 3> m_E_r;
  StorageT<ValueT, 3> m_E_z;
  StorageT<ValueT, 3> m_E_phi;

  ValueT m_t_min, m_t_max;
  ValueT m_r_min, m_r_max;
  ValueT m_z_min, m_z_max;
};

using DenseWeightingField = WeightingField<DenseNDArray>;

#endif
