#ifndef __DENSE_NDVECARRAY_HH
#define __DENSE_NDVECARRAY_HH

#include <iostream>

#include <memory>
#include <vector>
#include <span>

#include "Vector.hh"

template <class T, std::size_t dims, std::size_t vec_dims>
class DenseNDVecArray {

public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = std::span<T, vec_dims>;
  
private:
  using data_t = std::vector<T>;
  using stride_t = Vector<std::size_t, dims + 1>;

public:
  
  DenseNDVecArray(const shape_t& shape, const T& value) : m_offset(0), m_shape(shape) {
    m_strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), m_strides.begin() + 1, std::multiplies<std::size_t>());
    m_strides *= vec_dims;

    // reverse strides so that slowest-changing index is the last one
    std::reverse(std::execution::unseq, m_strides.begin(), m_strides.end());

    std::cout << "have the following final strides:" << std::endl;
    for(auto cur : m_strides) {
      std::cout << cur << std::endl;
    }

    std::cout << "volume = " << GetVolume() << std::endl;
    
    // reserve the required memory
    m_data -> resize(GetVolume(), value);
  }

  view_t operator[](const ind_t& ind) {
    std::size_t flat_ind = std::inner_product(ind.cbegin(), ind.cend(), m_strides.begin() + 1, m_offset);
    
    std::cout << "flat_ind = " << flat_ind << std::endl;
    
    return *m_data[flat_ind];
  }
  
  const std::size_t GetVolume() const {return m_strides[0];} 
    
private:
  
  std::shared_ptr<data_t> m_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  shape_t m_shape;
  
};

#endif
