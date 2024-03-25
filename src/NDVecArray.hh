#pragma once

#include <iostream> // for debug only---remove at the end

#include <fstream>
#include <memory>
#include <vector>
#include <span>

#include "Serialization.hh"
#include "Vector.hh"

// Forward declaration of (de)serializer
namespace stor {
  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct NDVecArrayStreamer;
}

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
class NDVecArray {

private:
  friend struct stor::NDVecArrayStreamer<T, dims, vec_dims>;
  
public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = VectorView<T, vec_dims>;
  
private:
  using data_t = std::vector<T>;
  using stride_t = Vector<std::size_t, dims + 1>;

  // constructor used by deserializer
  // TODO: to be deprecated and removed
  NDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, data_t&& data) :
    m_shape(shape), m_strides(strides), m_offset(offset), m_data(std::make_shared<data_t>(data)) { }

  // used in creation of view
  NDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, std::shared_ptr<data_t> data) :
    m_shape(shape), m_strides(strides), m_offset(offset), m_data(data) { }
  
public:
  
  NDVecArray(const shape_t& shape, const T& value) : m_offset(0), m_shape(shape) {
    m_strides = ComputeStrides(shape);
    
    // reserve the required memory
    m_data = std::make_shared<data_t>(GetVolume(), value);
  }   

  // element values are undefined after this operation (if size is increased), need to be set explicitly again
  void resize(const shape_t& new_shape) {
    m_shape = new_shape;
    m_strides = ComputeStrides(new_shape);
    m_offset = 0;
    m_data -> reserve(GetVolume());
  }
    
  // Single-element access
  view_t operator[](const ind_t& ind) {
    return view_t(m_data -> begin() + ComputeFlatInd(ind), vec_dims);
  }

  const view_t operator[](const ind_t& ind) const {
    return view_t(m_data -> begin() + ComputeFlatInd(ind), vec_dims);
  }
  
  // Array view access
  NDVecArray<T, dims, vec_dims> View(const ind_t& start_ind, const ind_t& end_ind) const {
    shape_t view_shape = end_ind - start_ind;
    stride_t view_strides = ComputeStrides(view_shape);        
    std::size_t view_offset = std::inner_product(start_ind.cbegin(), start_ind.cend(), m_strides.begin() + 1, m_offset);

    return NDVecArray<T, dims, vec_dims>(view_shape, view_strides, view_offset, m_data);
  }

  bool IsZero(const ind_t& ind) const {
    for(scalar_t& cur : this -> operator[](ind)) {
      if(cur != 0) {
	return false;
      }
    }
    return true;
  }
  
  const shape_t GetShape() const {return m_shape;}
  const std::size_t GetVolume() const {return m_strides[0];}
  const std::size_t GetNumberElements() const {return m_strides[0] / vec_dims;}
  
private:

  std::size_t ComputeFlatInd(const ind_t& ind) const {
    return std::inner_product(ind.cbegin(), ind.cend(), m_strides.begin() + 1, m_offset);
  }
  
  std::size_t ComputeFlatInd(const ind_t& ind) const requires(dims == 2) {
    return m_offset + ind[0] * m_strides[1] + ind[1] * m_strides[2];
  }

  std::size_t ComputeFlatInd(const ind_t& ind) const requires(dims == 3) {
    return m_offset + ind[0] * m_strides[1] + ind[1] * m_strides[2] + ind[2] * m_strides[3];
  }
  
  static stride_t ComputeStrides(const shape_t& shape) {
    stride_t strides;
    strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), strides.begin() + 1, std::multiplies<std::size_t>());
    strides *= vec_dims;

    // reverse strides so that slowest-changing index is the last one
    std::reverse(std::execution::unseq, strides.begin(), strides.end());
  
    return strides;
  }

  static std::size_t ComputeVolume(const shape_t& shape) {
    return ComputeStrides(shape)[0];
  }
  
private:
  
  std::shared_ptr<data_t> m_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  shape_t m_shape;
  
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArrayZeroAware : public NDVecArray<T, dims, vec_dims> {

  // with faster IsZero overload / bookkeeping of fraction of zero'ed elements -> to be used in compression step
  
};
