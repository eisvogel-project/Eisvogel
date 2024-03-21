#ifndef __DENSE_NDVECARRAY_HH
#define __DENSE_NDVECARRAY_HH

#include <iostream>

#include <memory>
#include <vector>
#include <span>

#include "Vector.hh"

template <class T, std::size_t vec_dims>
struct VectorView : public std::span<T, vec_dims> {

  template <class It>
  constexpr VectorView(It first, std::size_t count) : std::span<T, vec_dims>(first, count) { }; 

  VectorView& operator=(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.cbegin(), vec_dims, this -> begin());
    return *this;
  }
};

template <class T, std::size_t dims, std::size_t vec_dims>
class DenseNDVecArray {

public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = VectorView<T, vec_dims>;
  
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
    
    // reserve the required memory
    m_data = std::make_shared<data_t>(GetVolume(), value);
  }

  view_t operator[](const ind_t& ind) {
    std::size_t flat_ind = std::inner_product(ind.cbegin(), ind.cend(), m_strides.begin() + 1, m_offset);        
    return view_t(m_data -> begin() + flat_ind, vec_dims);
  }
  
  const std::size_t GetVolume() const {return m_strides[0];} 
    
private:
  
  std::shared_ptr<data_t> m_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  shape_t m_shape;
  
};

#endif
