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
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  class NDVecArrayStreamer;
}

template <typename T, std::size_t vec_dims>
struct VectorView : public std::span<T, vec_dims> {

  template <class It>
  constexpr VectorView(It first) : std::span<T, vec_dims>(first, vec_dims) { }; 

  VectorView& operator=(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.cbegin(), vec_dims, this -> begin());
    return *this;
  }
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArray {

  static_assert(dims > 0);
  static_assert(vec_dims > 0);
  
private:
  friend class stor::NDVecArrayStreamer<NDVecArray, T, dims, vec_dims>;
  
public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = VectorView<T, vec_dims>;
  
private:
  using data_t = std::vector<T>;
  using stride_t = Vector<std::size_t, dims>;

  // used in creation of view
  NDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, std::shared_ptr<data_t> data) :
    m_data(data), m_strides(strides), m_offset(offset), m_shape(shape) {
    m_number_elements =ComputeNumberElements(shape);
    m_volume = ComputeVolume(shape);
  }
  
public:
  
  NDVecArray(const shape_t& shape, const T& value) : m_offset(0), m_shape(shape) {
    m_strides = ComputeStrides(shape);
    m_number_elements = ComputeNumberElements(shape);
    m_volume = ComputeVolume(shape);
    
    // reserve the required memory
    m_data = std::make_shared<data_t>(GetVolume(), value);
  }   

  // element values are undefined after this operation (if size is increased), need to be set explicitly again
  void resize(const shape_t& new_shape, const T& value) {
    m_shape = new_shape;
    m_strides = ComputeStrides(new_shape);
    m_offset = 0;
    m_data -> resize(GetVolume(), value);
  }
    
  // Single-element access
  view_t operator[](const ind_t& ind) {
    return view_t(m_data -> begin() + ComputeFlatInd(ind));
  }

  const view_t operator[](const ind_t& ind) const {
    return view_t(m_data -> begin() + ComputeFlatInd(ind));
  }
  
  // Array view access
  NDVecArray<T, dims, vec_dims> View(const ind_t& start_ind, const ind_t& end_ind) const {
    std::size_t view_offset = ComputeFlatInd(start_ind);
    return NDVecArray<T, dims, vec_dims>(end_ind - start_ind, // shape of the view
					 m_strides, view_offset, m_data);
  }
  
  bool IsNull(const ind_t& ind) const {
    for(T& cur : this -> operator[](ind)) {
      if(cur != 0) {
	return false;
      }
    }
    return true;
  }
  
  const shape_t GetShape() const {return m_shape;}
  const std::size_t GetVolume() const {return m_volume;}
  const std::size_t GetNumberElements() const {return m_number_elements;}

  // determines whether array with shape `arr_shape` can be concatenated with an array with shape `other_shape` along `axis`
  static bool ShapeAllowsConcatenation(const shape_t& arr_shape, const shape_t& other_shape, std::size_t axis) {
    
    // array axis out of bounds
    if(axis >= dims) {
      return false;
    }

    // dimension along all directions need to match with the exception of the concatenation axis
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(ind == axis) {	continue; }
      if(arr_shape[ind] != other_shape[ind]) {
	return false;
      }
    }
    return true;
  }
  
  // -----------------------------
  // operations on `NDVecArrays`
  // -----------------------------

  // Copy range from other array
  // friend void CopyRange(const NDVecArray<T, dims, vec_dims>& source, const ind_t& source_start_ind, const ind_t& source_end_ind,
  //                       NDVecArray<T, dims, vec_dims>& target, const ind_t& target_start_ind) {  }
  
private:

  std::size_t ComputeFlatInd(const ind_t& ind) const {
    return std::inner_product(ind.cbegin(), ind.cend(), m_strides.cbegin(), m_offset);
  }

  std::size_t ComputeFlatInd(const ind_t& ind) const requires(dims == 1) {
    return m_offset + ind[0] * m_strides[0];
  }
  
  std::size_t ComputeFlatInd(const ind_t& ind) const requires(dims == 2) {
    return m_offset + ind[0] * m_strides[0] + ind[1] * m_strides[1];
  }

  std::size_t ComputeFlatInd(const ind_t& ind) const requires(dims == 3) {
    return m_offset + ind[0] * m_strides[0] + ind[1] * m_strides[1] + ind[2] * m_strides[2];
  }

  static stride_t ComputeStrides(const shape_t& shape) requires(dims == 1) {
    stride_t strides{vec_dims};
    return strides;
  }

  static stride_t ComputeStrides(const shape_t& shape) requires(dims == 2) {
    stride_t strides{shape[1] * vec_dims, vec_dims};
    return strides;
  }
  
  static stride_t ComputeStrides(const shape_t& shape) requires(dims > 2) {
    
    stride_t strides;
    std::size_t stride_accum = 1;
    for(std::size_t ind = dims - 1; ind != (std::size_t)(-1); ind--) {
      strides[ind] = stride_accum;
      stride_accum *= shape[ind];
    }    
    strides *= vec_dims;
    
    return strides;
  }

  static std::size_t ComputeNumberElements(const shape_t& shape) {
    return std::accumulate(shape.cbegin(), shape.cend(), 1, std::multiplies<std::size_t>());
  }

  static std::size_t ComputeVolume(const shape_t& shape) {
    return ComputeNumberElements(shape) * vec_dims;
  }
  
private:
  
  std::shared_ptr<data_t> m_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  std::size_t m_number_elements;
  std::size_t m_volume;
  shape_t m_shape;
  
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArrayNullAware : public NDVecArray<T, dims, vec_dims> {

  // with faster IsNull overload / bookkeeping of fraction of null'ed elements -> to be used in compression step
  
};