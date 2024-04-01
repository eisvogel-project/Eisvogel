#include <iostream>
#include <cassert>

// constructors
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const shape_t& shape) : m_offset(0), m_shape(shape), m_owns_data(true) {

  std::cout << "NDVecArray normal constructor uninitialized" << std::endl;
  
  m_strides = ComputeStrides(m_shape);
  UpdateShapeAttributes(m_shape);
  
  // reserve the required memory
  m_data = std::make_shared<data_t>(GetVolume());
}

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const shape_t& shape, const T& value) :
  NDVecArray<T, dims, vec_dims>(shape) {

  std::cout << "NDVecArray normal constructor initialized" << std::endl;
    
  // initialize properly
  this -> operator=(value);
}   

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const NDVecArray<T, dims, vec_dims>& other) :  
  m_offset(other.m_offset), m_shape(other.m_shape), m_strides(other.m_strides), m_number_elements(other.m_number_elements), m_volume(other.m_volume), m_owns_data(true) {

  std::cout << "NDVecArray copy constructor" << std::endl;
  
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

// copy-assignment operator
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>& NDVecArray<T, dims, vec_dims>::operator=(const NDVecArray<T, dims, vec_dims>& other) {
  
  resize(other.m_shape);

  // copy a region of memory that is guaranteed to be contiguous in both `other` and `this`
  auto copy_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) {
    
    std::size_t offset = ComputeFlatInd(chunk_begin);
    auto src = other.m_data -> begin() + offset;
    auto dest = m_data -> begin() + offset;
    
    std::copy_n(std::execution::unseq, src, m_shape[dims - 1] * vec_dims, dest);
  };

  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = m_shape[dims - 1];
  index_loop_over_array_chunks(*this, chunk_size, copy_chunk_contiguous);
  
  return *this;
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
NDVecArray<T, dims, vec_dims>& NDVecArray<T, dims, vec_dims>::operator=(const T& other) {

  // TODO: fix this, right now this would write across the end of a view
  std::fill(m_data -> begin(), m_data -> end(), other);
  return *this;
}

// loop over elements
template <typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void NDVecArray<T, dims, vec_dims>::loop_over_elements(CallableT&& worker) const {

  // manual handling of the loop over the contiguous memory region
  auto loop_over_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) {

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
  index_loop_over_array_chunks(*this, chunk_size, loop_over_chunk_contiguous);    
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
