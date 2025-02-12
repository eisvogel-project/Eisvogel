#pragma once

#include <memory>
#include <limits>
#include "NDVecArray.hh"
#include "MemoryUtils.hh"

namespace stor {

  enum class StreamerMode : std::size_t {
    automatic = 0,
    dense = 1,
    null_suppressed = 2
  };

  // Controls whether or not the on-disk representation of the array can be modified or not
  // Modifications may be disallowed for various reasons, e.g. because the on-disk array has been marked as `final`,
  // or because unequal serialization chunks have been used
  enum class AccessMode : std::size_t {
    modification_allowed = 0,
    modification_not_allowed = 1
  };
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  struct NDVecArrayStreamerMetadata;

  static constexpr std::size_t INFTY = std::numeric_limits<std::size_t>::max();
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  class NDVecArrayStreamer {

  private:
    friend struct NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>;

    using type = ArrayT<T, dims, vec_dims>;
    using data_t = typename type::data_t;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;

    // The type of the serialized data
    using ser_type = Traits<T>::ser_type;

    // Note: if the calls to std::vector<ser_type>::resize ever become the bottleneck due to
    // repeated zero initialization of all ements, replace with something else
    // (the initial values are overwritten anyways on any read from disk, so no
    // explicit zero initialization is ever required)
    using buffer_t = std::vector<ser_type, no_init_alloc<ser_type>>;
    
  public:

    NDVecArrayStreamer(std::size_t initial_buffer_size = 10000);

    void mark_as_final(std::fstream& stream);
    
    void serialize(std::fstream& stream, const type& val, const shape_t& chunk_size, const StreamerMode& mode);
    void deserialize(std::fstream& stream, type& val);

    void append_slice(std::fstream& stream, const type& chunk, std::size_t axis, const StreamerMode& mode);

  private:

    void serialize_all_chunks(std::fstream& stream, const type& val, const shape_t& chunk_size, const StreamerMode& mode);
    void serialize_all_chunks_dense(std::fstream& stream, const type& val, const shape_t& chunk_size);
    void serialize_all_chunks_null_suppressed(std::fstream& stream, const type& val, const shape_t& chunk_size);

  private:

    void write_buffer(std::size_t num_elems, std::fstream& stream);
    void read_into_buffer(std::size_t num_elems, std::fstream& stream);
    
  private:

    std::shared_ptr<buffer_t> m_ser_buffer;

  };
}

#include "NDVecArrayStreamer.hxx"
