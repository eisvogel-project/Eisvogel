#ifndef __DENSE_NDVECARRAY_HH
#define __DENSE_NDVECARRAY_HH

#include <iostream>

#include <memory>
#include <vector>
#include <span>

#include "Serialization.hh"
#include "Vector.hh"

template <typename T, std::size_t vec_dims>
struct VectorView : public std::span<T, vec_dims> {

  template <class It>
  constexpr VectorView(It first, std::size_t count) : std::span<T, vec_dims>(first, count) { }; 

  VectorView& operator=(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.cbegin(), vec_dims, this -> begin());
    return *this;
  }
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class DenseNDVecArray {

public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = VectorView<T, vec_dims>;
  
private:
  using data_t = std::vector<T>;
  using stride_t = Vector<std::size_t, dims + 1>;

private:
  friend struct stor::Traits<DenseNDVecArray<T, dims, vec_dims>>;

  // constructor used by deserializer
  DenseNDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, data_t&& data) :
    m_shape(shape), m_strides(strides), m_offset(offset), m_data(std::make_shared<data_t>(data)) { }

  DenseNDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, std::shared_ptr<data_t> data) :
    m_shape(shape), m_strides(strides), m_offset(offset), m_data(data) { }
  
public:
  
  DenseNDVecArray(const shape_t& shape, const T& value) : m_offset(0), m_shape(shape) {
    m_strides = ComputeStrides(shape);
    
    // reserve the required memory
    m_data = std::make_shared<data_t>(GetVolume(), value);
  }   
  
  // Single-element access
  view_t operator[](const ind_t& ind) {
    std::size_t flat_ind = std::inner_product(ind.cbegin(), ind.cend(), m_strides.begin() + 1, m_offset);
    return view_t(m_data -> begin() + flat_ind, vec_dims);
  }

  view_t operator[](const ind_t& ind) requires(dims == 2) {
    std::size_t flat_ind = m_offset + ind[0] * m_strides[1] + ind[1] * m_strides[2];
    return view_t(m_data -> begin() + flat_ind, vec_dims);
  }
  
  view_t operator[](const ind_t& ind) requires(dims == 3) {
    std::size_t flat_ind = m_offset + ind[0] * m_strides[1] + ind[1] * m_strides[2] + ind[2] * m_strides[3];
    return view_t(m_data -> begin() + flat_ind, vec_dims);
  }

  // Array view access
  DenseNDVecArray<T, dims, vec_dims> View(const ind_t& start_ind, const ind_t& end_ind) {
    shape_t view_shape = end_ind - start_ind;
    stride_t view_strides = ComputeStrides(view_shape);        
    std::size_t view_offset = std::inner_product(start_ind.cbegin(), start_ind.cend(), m_strides.begin() + 1, m_offset);

    return DenseNDVecArray<T, dims, vec_dims>(view_shape, view_strides, view_offset, m_data);
  }
  
  const std::size_t GetVolume() const {return m_strides[0];} 

private:

  static stride_t ComputeStrides(const shape_t& shape) {
    stride_t strides;
    strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), strides.begin() + 1, std::multiplies<std::size_t>());
    strides *= vec_dims;

    // reverse strides so that slowest-changing index is the last one
    std::reverse(std::execution::unseq, strides.begin(), strides.end());
  
    return strides;
  }
  
private:
  
  std::shared_ptr<data_t> m_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  shape_t m_shape;
  
};

namespace stor {

  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct Traits<DenseNDVecArray<T, dims, vec_dims>> {
    using type = DenseNDVecArray<T, dims, vec_dims>;
    using data_t = typename type::data_t;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<shape_t>::serialize(stream, val.m_shape);
      Traits<stride_t>::serialize(stream, val.m_strides);
      Traits<std::size_t>::serialize(stream, val.m_offset);
      Traits<data_t>::serialize(stream, *val.m_data);
    }

    static type deserialize(std::iostream& stream) {
      shape_t shape = Traits<shape_t>::deserialize(stream);
      stride_t strides = Traits<stride_t>::deserialize(stream);
      std::size_t offset = Traits<std::size_t>::deserialize(stream);
      data_t data = Traits<data_t>::deserialize(stream);

      return type(shape, strides, offset, std::move(data));
    }
  };  
}

#endif
