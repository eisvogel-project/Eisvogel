#include "DenseWeightingFieldProvider.hh"

DenseWeightingFieldProvider::DenseWeightingFieldProvider(){
  
}

DenseWeightingFieldProvider::value_t DenseWeightingFieldProvider::get_E_r(value_t t, value_t r, value_t z) {
  return 1.0;
}

value_t DenseWeightingFieldProvider::get_E_r(index_t ind_t, index_t ind_r, index_t ind_z) {
  return 1.0;
}

value_t DenseWeightingFieldProvider::get_E_z(value_t t, value_t r, value_t z) {
  return 2.0;
}

value_t DenseWeightingFieldProvider::get_E_z(index_t ind_t, index_t ind_r, index_t ind_z) {
  return 2.0;
}

value_t DenseWeightingFieldProvider::get_E_phi(value_t t, value_t r, value_t z) {
  return 3.0;
}

value_t DenseWeightingFieldProvider::get_E_phi(index_t ind_t, index_t ind_r, index_t ind_z) {
  return 3.0;
}
