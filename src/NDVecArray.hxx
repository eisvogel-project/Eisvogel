#include <iostream>
#include <cassert>

// constructors
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const shape_t& shape) : m_owns_data(true), m_offset(0), m_shape(shape) {
  
  m_strides = ComputeStrides(m_shape);
  UpdateShapeAttributes(m_shape);
  
  // reserve the required memory
  m_data = std::make_shared<data_t>(GetVolume());
}

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const shape_t& shape, const T& value) :
  NDVecArray<T, dims, vec_dims>(shape) {
    
  // initialize properly
  this -> operator=(value);
}   

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(std::size_t len) requires(dims == 1) : NDVecArray(shape_t{len}) { };

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(std::size_t len, const T& value) requires(dims == 1) : NDVecArray(shape_t{len}, value) { };

// coppy constructor
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const NDVecArray<T, dims, vec_dims>& other) :  
  m_owns_data(true), m_strides(other.m_strides), m_offset(other.m_offset), m_number_elements(other.m_number_elements), m_volume(other.m_volume), m_shape(other.m_shape)  {
  
  // reserve the required memory
  m_data = std::make_shared<data_t>(GetVolume());
  
  // copy data
  this -> operator=(other);
}

// resizing operations
template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::resize(const shape_t& new_shape) {

  // Need to own our own data for this to make sense
  assert(m_owns_data);
  
  m_shape = new_shape;
  m_strides = ComputeStrides(new_shape);
  m_offset = 0;

  UpdateShapeAttributes(m_shape);
  
  m_data -> resize(GetVolume());
}

template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::resize(const shape_t& new_shape, const T& value) {  
  resize(new_shape);
  this -> operator=(value);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::resize(std::size_t new_len) requires(dims == 1) {
  resize(shape_t{new_len});
}

template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::clear() {
  assert(m_owns_data);
  std::fill(std::execution::unseq, m_data -> begin(), m_data -> end(), (T)0);
}

// copy-assignment operator
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>& NDVecArray<T, dims, vec_dims>::operator=(const NDVecArray<T, dims, vec_dims>& other) {
  
  resize(other.m_shape);

  // copy a region of memory that is guaranteed to be contiguous in both `other` and `this`
  auto copy_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {
    
    std::size_t offset = ComputeFlatInd(chunk_begin);
    auto src = other.m_data -> begin() + offset;
    auto dest = m_data -> begin() + offset;
    
    std::copy_n(std::execution::unseq, src, m_shape[dims - 1] * vec_dims, dest);
  };

  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = m_shape[dims - 1];
  IteratorUtils::index_loop_over_array_chunks(*this, chunk_size, copy_chunk_contiguous);
  
  return *this;
}

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>& NDVecArray<T, dims, vec_dims>::operator=(const T& other) {

  // fill a contiguous chunk of memory with the target value
  auto fill_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {
    
    std::size_t offset = ComputeFlatInd(chunk_begin);
    auto dest = m_data -> begin() + offset;
    
    std::fill_n(std::execution::unseq, dest, m_shape[dims - 1] * vec_dims, other);
  };

  Vector<std::size_t, dims> chunk_size(1u);
  chunk_size[dims - 1] = m_shape[dims - 1];
  IteratorUtils::index_loop_over_array_chunks(*this, chunk_size, fill_chunk_contiguous);

  return *this;
}

template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::fill_from(const NDVecArray<T, dims, vec_dims>& other,
					      const ind_t& input_start, const ind_t& input_end, const ind_t& output_start) {

  // make sure the full range is available in the target and source arrays
  assert(other.has_index(input_start));
  assert(other.has_index(input_end - 1));  // the end index is always exclusive
  assert(has_index(output_start));
  assert(has_index(output_start + (input_end - input_start) - 1));  // the end index is always exclusive

  // copy a region of memory that is guaranteed to be contiguous in both `other` and `this`
  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = input_end[dims - 1] - input_start[dims - 1];

  auto copy_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {
    
    auto src = other.m_data -> begin() + other.ComputeFlatInd(chunk_begin);
    auto dest = m_data -> begin() + ComputeFlatInd(output_start + (chunk_begin - input_start));
    
    std::copy_n(std::execution::unseq, src, chunk_size[dims - 1] * vec_dims, dest);
  };

  // iterate over chunks in the source array
  IteratorUtils::index_loop_over_chunks(input_start, input_end, chunk_size, copy_chunk_contiguous);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
bool NDVecArray<T, dims, vec_dims>::has_index(const ind_t& ind, const shape_t& shape) {
  for(std::size_t i = 0; i < dims; i++) {
    if(ind[i] >= shape[i]) {
      return false;
    }
  }
  return true;
}

template <typename T, std::size_t dims, std::size_t vec_dims>
bool NDVecArray<T, dims, vec_dims>::has_index(const ind_t& ind) const {
  return has_index(ind, m_shape);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void NDVecArray<T, dims, vec_dims>::Append(const NDVecArray<T, dims, vec_dims>& other) requires(axis == 0) {
  
  // Need to own our own data for this to make sense
  assert(m_owns_data);
  
  // some sanity checks
  static_assert(axis < dims);
  assert(ShapeAllowsConcatenation(m_shape, other.m_shape, axis));
  
  // new elements will be appended at the very back
  std::size_t app_offset = m_data -> size();

  // compute new shape and resize
  shape_t new_shape = m_shape;
  new_shape[axis] += other.m_shape[axis];
  resize(new_shape);
  
  // copy the new data
  auto it_app = m_data -> begin() + app_offset;
  std::copy(std::execution::unseq, other.m_data -> begin(), other.m_data -> end(), it_app);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void NDVecArray<T, dims, vec_dims>::Append(const NDVecArray<T, dims, vec_dims>& other) {
  static_assert(axis < dims);
  assert(m_owns_data);
  throw std::logic_error("Not implemented yet!");
}

template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::Append(const NDVecArray<T, dims, vec_dims>& other, std::size_t axis) {
  if(axis == 0) {
    Append<0>(other);
  }
  else {
    assert(axis < dims);
    assert(m_owns_data);
    throw std::logic_error("Not implemented yet!");
  }
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis_1, std::size_t axis_2>
void NDVecArray<T, dims, vec_dims>::SwapAxes(NDVecArray<T, dims, vec_dims>& dest) const {

  // make sure we have the correct shape
  shape_t swapped_shape = m_shape;
  std::swap(swapped_shape[axis_1], swapped_shape[axis_2]);
  dest.resize(swapped_shape);
  
  auto copy_swapped = [&](const Vector<std::size_t, dims>& ind, const view_t& element) {

    // prepare the index with the requested dimensions swapped
    ind_t swapped_ind = ind;
    std::swap(swapped_ind[axis_1], swapped_ind[axis_2]);
    
    dest[swapped_ind] = element;
  };
  index_loop_over_elements(copy_swapped);
}

// loop over elements
template <typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void NDVecArray<T, dims, vec_dims>::loop_over_elements(CallableT&& worker) const {  
  
  // manual handling of the loop over the contiguous memory region
  auto loop_over_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {

    // iterator to the beginning of this contiguous memory region
    auto it = m_data -> begin() + ComputeFlatInd(chunk_begin);

    // end of this contiguous memory region
    auto it_end = it + m_shape[dims - 1] * vec_dims;

    // call worker function on every element
    while(it != it_end) {
      worker(view_t(it));
      it += vec_dims;
    }    
  };

  // iterate over chunks that are contiguous in memory, i.e. one stride of the fastest-varying index
  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = m_shape[dims - 1];
  IteratorUtils::index_loop_over_array_chunks(*this, chunk_size, loop_over_chunk_contiguous);    
}

// Loop over (index, element) pairs
template <typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void NDVecArray<T, dims, vec_dims>::index_loop_over_elements(const ind_t& start_ind, const ind_t& end_ind, CallableT&& worker) const {

  Vector<std::size_t, dims> ind;

  // manual handling of the iteration over chunks that are contiguous in memory, i.e. one stride of the fastest-varying index
  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = end_ind[dims - 1] - start_ind[dims - 1];
  
  auto loop_over_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {
    
    // iterator to the beginning of this contiguous memory region
    auto it = m_data -> begin() + ComputeFlatInd(chunk_begin);

    // end of this contiguous memory region
    auto it_end = it + chunk_size[dims - 1] * vec_dims;

    // call worker function on every (index, element) pair
    ind = chunk_begin;
    while(it != it_end) {
      worker(ind, view_t(it));
      it += vec_dims;
      ind[dims - 1] += 1;
    }    
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_size, loop_over_chunk_contiguous);  
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void NDVecArray<T, dims, vec_dims>::index_loop_over_elements(CallableT&& worker) const {
  Vector<std::size_t, dims> start_ind(0);
  index_loop_over_elements(start_ind, m_shape, worker);
}

// Sequential access
// Note: this does not seem significantly faster than manual sequential access through operator[]
template <typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void NDVecArray<T, dims, vec_dims>::apply_sequential(const ind_t& outer_ind,
							       const std::vector<std::size_t>& inner_ind,
							       CallableT&& worker) {
  
  auto it_outer = m_data -> begin() + ComputeFlatInd(outer_ind);
  for(const std::size_t& cur_inner_ind : inner_ind) {
    worker(view_t(it_outer + cur_inner_ind * vec_dims));
  }
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
bool NDVecArray<T, dims, vec_dims>::ShapeAllowsConcatenation(const shape_t& arr_shape, const shape_t& other_shape) {
  return ShapeAllowsConcatenation(arr_shape, other_shape, axis);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
bool NDVecArray<T, dims, vec_dims>::ShapeAllowsConcatenation(const shape_t& arr_shape, const shape_t& other_shape, std::size_t axis) {
    
  // array axis out of bounds
  assert(axis < dims);
  
  // dimension along all directions need to match with the exception of the concatenation axis
  for(std::size_t ind = 0; ind < dims; ind++) {
    if(ind == axis) { continue; }
    if(arr_shape[ind] != other_shape[ind]) {
      return false;
    }
  }
  return true;
}

namespace stor {

  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct Traits<NDVecArray<T, dims, vec_dims>> {
    using type = NDVecArray<T, dims, vec_dims>;

    static void serialize_to_numpy(std::iostream& stream, const type& val) {
      
    }
  };  
}
