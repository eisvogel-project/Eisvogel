#ifndef __WEIGHTING_FIELD_PROVIDER_HH
#define __WEIGHTING_FIELD_PROVIDER_HH

class WeightingFieldProvider {

  using index_t = unsigned int;
  using value_t = double;

public:
 
  virtual value_t get_E_r(value_t t, value_t r, value_t z) = 0;
  virtual value_t get_E_z(value_t t, value_t r, value_t z) = 0;
  virtual value_t get_E_phi(value_t t, value_t r, value_t z) = 0;

  virtual value_t get_E_r(index_t ind_t, index_t ind_r, index_t ind_z) = 0;
  virtual value_t get_E_z(index_t ind_t, index_t ind_r, index_t ind_z) = 0;
  virtual value_t get_E_phi(index_t ind_t, index_t ind_r, index_t ind_z) = 0;

private:

  value_t rstart, rend;
  value_t zstart, zend;
  value_t tstart, tend;

};

#endif
