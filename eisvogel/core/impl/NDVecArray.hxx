#include <iostream>
#include <cassert>
#include <format>

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

  // std::cout << "in NDVecArray copy assignment operator" << std::endl;
  
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
void NDVecArray<T, dims, vec_dims>::fill(const ind_t& start, const ind_t& end, const T& value) {

  // make sure the full range to be filled is available in this array
  assert(has_index(start));
  assert(has_index(end));

  // fill the range from `start` to `end` in chunks that are guaranteed to be contiguous in memory
  Vector<std::size_t, dims> chunk_size(1);
  chunk_size[dims - 1] = end[dims - 1] - start[dims - 1];

  auto fill_chunk_contiguous = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>&) {
    auto dest = m_data -> begin() + ComputeFlatInd(chunk_begin);
    std::fill_n(std::execution::unseq, dest, chunk_size[dims - 1] * vec_dims, value);
  };

  // iterate over all such contiguous chunks
  IteratorUtils::index_loop_over_chunks(start, end, chunk_size, fill_chunk_contiguous);
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
template <std::size_t axis>
void NDVecArray<T, dims, vec_dims>::fill_from(const NDVecArray<T, 1, vec_dims>& other,
					      const ind_t& output_start) {
  static_assert(axis < dims);
  assert(m_owns_data);
  throw std::logic_error("Not implemented yet!");
}

template <typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void NDVecArray<T, dims, vec_dims>::fill_from(const NDVecArray<T, 1, vec_dims>& other,
					      const ind_t& output_start) requires(axis == dims - 1) {
  assert(m_owns_data);
  assert(other.owns_data());

  // The `other` array is 1-dimensional
  std::size_t other_shape = other.GetShape()[0];
  
  // make sure the full range is available
  ind_t output_end = output_start;
  output_end[axis] += other_shape - 1;  // the upper range is exclusive, as always

  // std::cout << "output_start = " << output_start << std::endl;
  // std::cout << "output_end = " << output_end << std::endl;
  // std::cout << "own shape = " << GetShape() << std::endl;
  // std::cout << "other shape = " << other.GetShape() << std::endl;
  
  assert(has_index(output_start));
  assert(has_index(output_end));

  // `other` fits into this array as a single contiguous block, can simply copy
  auto src = other.data();
  auto dest = m_data -> begin() + ComputeFlatInd(output_start);
  std::copy_n(std::execution::unseq, src, vec_dims * other_shape, dest);
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

// Serialization to numpy binary file
namespace stor {

  namespace Numpy {

    template <typename T>
    struct Traits;

    template <>
    struct Traits<float> {
      // Everything is serialized in network byte-order, which is defined to be big-endian
      static constexpr std::string dtype = ">f4";
      using ser_type = uint32_t;
    };
  }
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct Traits<NDVecArray<T, dims, vec_dims>> {
    using type = NDVecArray<T, dims, vec_dims>;
    using ser_type = typename Numpy::Traits<T>::ser_type;

    // See description of Numpy file format (v1.0) at https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html
    static void serialize_to_numpy(std::iostream& stream, const type& val) {

      assert(val.m_owns_data);
      
      // Magic string to signify a numpy file
      const std::string init = "\x93NUMPY";
      serialize_string(stream, init);

      // File format version
      const char major_version = 0x01;
      const char minor_version = 0x00;
      stream.write(&major_version, 1);
      stream.write(&minor_version, 1);

      // Header with some metadata
      std::string header = std::format("{{'descr': '{}', 'fortran_order': False, 'shape': {}, }}",
				       Numpy::Traits<T>::dtype, shape_string(val));

      std::size_t header_length = header.size();
      
      // Now need to determine how many "\x20" to append to the header string to make the total file header
      // (from the beginning) to be divisible by 64      
      std::size_t length_before_padding = init.size() + sizeof(major_version) + sizeof(minor_version) + 2 + header_length;
      std::size_t padding_length = (length_before_padding / 64 + 1) * 64 - length_before_padding;
      header.append(padding_length - 1, 0x20);
      header.append("\n");

      std::size_t padded_header_length = header.size();
      
      // Serialize header length and header data
      serialize_short_little_endian(stream, padded_header_length);
      serialize_string(stream, header);

      // Serialize array data
      serialize_data(stream, val);
    }

  private:
    
    static void serialize_string(std::iostream& stream, const std::string& str) {
      stream.write(str.data(), str.size());
    }

    static void serialize_short_little_endian(std::iostream& stream, unsigned short int val) {
      char lo = (val >> 0) & 0xFF;
      char hi = (val >> 8) & 0xFF;      
      stream.write(&lo, sizeof(lo));
      stream.write(&hi, sizeof(hi));      
    }
    
    static std::string shape_string(const type& val) {
      std::string retval = "(";
      for(std::size_t i = 0; i < dims; i++) {
	retval += std::to_string(val.m_shape[i]) + ", ";
      }
      retval += std::to_string(vec_dims) + ")";
      return retval;
    }

    static void serialize_data(std::iostream& stream, const type& val, std::size_t block_size = 10000) {
      std::size_t vec_size = val.m_data -> size();
      std::size_t vec_ind = 0;

      while(vec_ind < vec_size) {
	std::vector<ser_type> outbuf(std::min(vec_size - vec_ind, block_size));
	for(std::size_t buf_ind = 0; (buf_ind < block_size) && (vec_ind < vec_size); buf_ind++, vec_ind++) {
	  ser_type ser_val = reinterpret_cast<const ser_type&>((*val.m_data)[vec_ind]);
	  outbuf[buf_ind] = htonl(ser_val);
	}
	stream.write((char*)outbuf.data(), sizeof(outbuf[0]) * outbuf.size());
      }
    }
  };  
}
