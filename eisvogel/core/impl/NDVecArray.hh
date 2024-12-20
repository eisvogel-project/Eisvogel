#pragma once

#include <iostream> // for debug only---remove at the end

#include <cassert>
#include <fstream>
#include <memory>
#include <vector>

#include "IteratorUtils.hh"
#include "Vector.hh"
#include "VectorView.hh"
#include "MemoryUtils.hh"

// Forward declaration of (de)serializer
namespace stor {
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  class NDVecArrayStreamer;
}

// Forward declaration of interpolation routine for NDVecArrays
template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArray;

namespace Interpolation::Kernel {
  struct Keys;
  struct Linear;
}

namespace Interpolation {
  
  template <class KernelT, typename T, std::size_t dims, std::size_t vec_dims>
  void interpolate(const NDVecArray<T, dims, vec_dims>& arr, NDVecArray<T, 1, vec_dims>& retval,
		   const Vector<scalar_t, dims-1>& outer_inds,
		   scalar_t inner_ind_start, scalar_t inner_ind_delta, std::size_t num_samples);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArray {

  static_assert(dims > 0);
  static_assert(vec_dims > 0);

private:
  friend struct stor::Traits<NDVecArray<T, dims, vec_dims>>;
  
private:

  // Streamer becomes a friend for (de)serialization
  friend class stor::NDVecArrayStreamer<NDVecArray, T, dims, vec_dims>;

  // Interpolators become friends for efficiency
  // (Note: need to manually list kernels because ISO C++20 does not allow partial template specializations as friends)
  friend void
  Interpolation::interpolate<Interpolation::Kernel::Keys, T, dims, vec_dims>(const NDVecArray<T, dims, vec_dims>&,
									     NDVecArray<T, 1, vec_dims>&,
									     const Vector<scalar_t, dims-1>&,
									     scalar_t, scalar_t, std::size_t);

  friend void
  Interpolation::interpolate<Interpolation::Kernel::Linear, T, dims, vec_dims>(const NDVecArray<T, dims, vec_dims>&,
									       NDVecArray<T, 1, vec_dims>&,
									       const Vector<scalar_t, dims-1>&,
									       scalar_t, scalar_t, std::size_t);
  
public:
  using ind_t = Vector<std::size_t, dims>;
  using shape_t = Vector<std::size_t, dims>;
  using view_t = VectorView<T, vec_dims>;
  
private:
  using data_t = std::vector<T, no_init_alloc<T>>;
  using stride_t = Vector<std::size_t, dims>;

  // used in creation of view
  NDVecArray(const shape_t& shape, const stride_t& strides, const std::size_t offset, std::shared_ptr<data_t> data) :
    m_data(data), m_owns_data(false), m_strides(strides), m_offset(offset),
    m_number_elements(ComputeNumberElements(shape)), m_volume(ComputeVolume(shape)), m_shape(shape) { }
  
public:

  // Creates an array with the specified `shape` that gets initialized with `value`
  NDVecArray(const shape_t& shape, const T& value);
  NDVecArray(std::size_t len, const T& value) requires(dims == 1);

  // Creates an array with the specified `shape` that remains uninitialized
  NDVecArray(const shape_t& shape);
  NDVecArray(std::size_t len) requires(dims == 1);

  // copy constructor
  NDVecArray(const NDVecArray<T, dims, vec_dims>& other);
  
  // copy-assignment operators (these copy the full array and thus also change the shape of the destination)
  // requires that the target `NDVecArray` owns its data (i.e. is not a view of another array)
  NDVecArray<T, dims, vec_dims>& operator=(const NDVecArray<T, dims, vec_dims>& other);

  // fills this array with the value `other`; also works if this array actually is only a view and doesn't own its data
  NDVecArray<T, dims, vec_dims>& operator=(const T& other);

  // fill the region from `start` to `end` in this array with `value`
  void fill(const ind_t& start, const ind_t& end, const T& value);
  
  // copy and fill part of the full array from `other` without changing the shape
  void fill_from(const NDVecArray<T, dims, vec_dims>& other,
		 const ind_t& input_start, const ind_t& input_end,
		 const ind_t& output_start);

  // Fill a 1-dimensional column along the direction of `axis`
  template <std::size_t axis>
  void fill_from(const NDVecArray<T, 1, vec_dims>& other,
		 const ind_t& output_start);

  template <std::size_t axis>
  void fill_from(const NDVecArray<T, 1, vec_dims>& other,
		 const ind_t& output_start) requires(axis == dims - 1);
  
  // element values are undefined after this operation (if size is increased), need to be set explicitly again
  void resize(const shape_t& new_shape);
  void resize(std::size_t new_len) requires(dims == 1);
  void resize(const shape_t& new_shape, const T& value);

  // sets all entries of this array to zero; requires that data is owned
  void clear();

  bool index_within_bounds(const ind_t& ind) const {
    for(std::size_t i = 0; i < dims; i++) {
      if(ind[i] >= m_shape[i]) {
	return false;
      }
    }
    return true;
  }
  
  // Single-element access, no bounds checking
  view_t operator[](const ind_t& ind) {
    assert(index_within_bounds(ind));
    return view_t(m_data -> begin() + ComputeFlatInd(ind));
  }

  const view_t operator[](const ind_t& ind) const {
    assert(index_within_bounds(ind));
    return view_t(m_data -> begin() + ComputeFlatInd(ind));
  }

  view_t operator[](const std::size_t ind) requires(dims == 1) {
    assert(index_within_bounds(ind));
    return view_t(m_data -> begin() + vec_dims * ind);
  }

  // T& operator[](const ind_t& ind) requires(vec_dims == 1) {
  //   assert(index_within_bounds(ind));
  //   return *(m_data -> begin() + ComputeFlatInd(ind));
  // }

  // T& operator[](const std::size_t ind) requires((dims == 1) && (vec_dims == 1)) {
  //   assert(index_within_bounds(ind));
  //   return m_data[ind];
  // }
  
  // Sequential access
  template <class CallableT>
  constexpr void apply_sequential(const ind_t& outer_ind, const std::vector<std::size_t>& inner_ind,
				  CallableT&& worker);
  
  // Loop over elements
  template <class CallableT>
  constexpr void loop_over_elements(CallableT&& worker) const;

  // Loop over (index, element) pairs
  template <class CallableT>
  constexpr void index_loop_over_elements(CallableT&& worker) const;

  template <class CallableT>
  constexpr void index_loop_over_elements(const ind_t& start_ind, const ind_t& end_ind, CallableT&& worker) const;
  
  // Array view access
  NDVecArray<T, dims, vec_dims> View(const ind_t& start_ind, const ind_t& end_ind) const {
    std::size_t view_offset = ComputeFlatInd(start_ind);
    return NDVecArray<T, dims, vec_dims>(end_ind - start_ind, // shape of the view
					 m_strides, view_offset, m_data);
  }
  
  bool IsNull(const ind_t& ind) const;
  bool IsAllNull() const;
  
  const shape_t& GetShape() const {return m_shape;}
  std::size_t GetVolume() const {return m_volume;}
  std::size_t GetNumberElements() const {return m_number_elements;}

  // determines whether array with shape `arr_shape` can be concatenated with an array with shape `other_shape` along `axis`
  template <std::size_t axis>
  static bool ShapeAllowsConcatenation(const shape_t& arr_shape, const shape_t& other_shape);
  
  static bool ShapeAllowsConcatenation(const shape_t& arr_shape, const shape_t& other_shape, std::size_t axis);

  // For appending another NDVecArray to this one
  template <std::size_t axis>
  void Append(const NDVecArray<T, dims, vec_dims>& other) requires(axis == 0);
  
  template <std::size_t axis>
  void Append(const NDVecArray<T, dims, vec_dims>& other);

  void Append(const NDVecArray<T, dims, vec_dims>& other, std::size_t axis);

  // Builds a version of this array at `dest` with `axis_1` and `axis_2` swapped
  template <std::size_t axis_1, std::size_t axis_2>
  void SwapAxes(NDVecArray<T, dims, vec_dims>& dest) const;
  
  // Other utilities
  bool has_index(const ind_t& ind) const;
  static bool has_index(const ind_t& ind, const shape_t& shape);

  bool owns_data() const {return m_owns_data;}
  auto data() const {return m_data -> begin();}
  
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

  static stride_t ComputeStrides(const shape_t&) requires(dims == 1) {
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

  void UpdateShapeAttributes(const shape_t& shape) {
    m_number_elements = ComputeNumberElements(shape);
    m_volume = ComputeVolume(shape);
  }
  
  static std::size_t ComputeNumberElements(const shape_t& shape) {
    return std::accumulate(shape.cbegin(), shape.cend(), (std::size_t)1u, std::multiplies<std::size_t>());
  }

  static std::size_t ComputeVolume(const shape_t& shape) {
    return ComputeNumberElements(shape) * vec_dims;
  }
  
private:
  
  std::shared_ptr<data_t> m_data;
  const bool m_owns_data;
  
  stride_t m_strides;
  std::size_t m_offset;
  std::size_t m_number_elements;
  std::size_t m_volume;
  shape_t m_shape;
  
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class NDVecArrayNullAware : public NDVecArray<T, dims, vec_dims> {

  // with faster IsNull and IsAllNull overload / bookkeeping of fraction of null'ed elements -> to be used in compression step
  
};

#include "NDVecArray.hxx"
