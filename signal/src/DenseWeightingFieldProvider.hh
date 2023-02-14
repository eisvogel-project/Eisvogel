#ifndef __DENSE_WEIGHTING_FIELD_PROVIDER_HH
#define __DENSE_WEIGHTING_FIELD_PROVIDER_HH

#include <unordered_map>
#include <tuple>
#include "WeightingFieldProvider.hh"

class DenseWeightingFieldProvider : public WeightingFieldProvider {

private:

  // Components indexed as (t, r, z)
  using index3d_t = std::tuple<index_t, index_t, index_t>;
  using storage_t = std::unordered_map<index3d_t, value_t>;

  storage_t E_r;
  storage_t E_z;
  storage_t E_phi;
};

#endif
