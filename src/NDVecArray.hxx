#include <iostream>

// constructors
template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const shape_t& shape, const T& value) : m_offset(0), m_shape(shape) {  
  m_strides = ComputeStrides(shape);
  m_number_elements = ComputeNumberElements(shape);
  m_volume = ComputeVolume(shape);
  
  // reserve the required memory
  m_data = std::make_shared<data_t>(GetVolume());
  
  // initialize properly
  this -> operator=(value);
}   

template <typename T, std::size_t dims, std::size_t vec_dims>
NDVecArray<T, dims, vec_dims>::NDVecArray(const NDVecArray<T, dims, vec_dims>& other) :  
  m_offset(other.m_offset), m_shape(other.m_shape), m_strides(other.m_strides), m_number_elements(other.m_number_elements), m_volume(other.m_volume) {
  
  // reserve the required memory
  m_data = std::make_shared<data_t>(GetVolume());
  
  // copy data
  this -> operator=(other);
}

// resizing operations
template <typename T, std::size_t dims, std::size_t vec_dims>
void NDVecArray<T, dims, vec_dims>::resize(const shape_t& new_shape) {
  m_shape = new_shape;
  m_strides = ComputeStrides(new_shape);
  m_offset = 0;
  
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

  std::cout << "NDVecArray copy-assignment" << std::endl;
  
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
NDVecArray<T, dims, vec_dims>& NDVecArray<T, dims, vec_dims>::operator=(const T& other) {
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
