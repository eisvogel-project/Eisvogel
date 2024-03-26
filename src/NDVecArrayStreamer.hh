#pragma once

#include <memory>
#include "NDVecArray.hh"

namespace stor {

  enum class StreamerMode : std::size_t {
    automatic = 0,
    dense = 1,
    zero_suppressed = 2
  };

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  struct NDVecArrayStreamerMetadata;
  
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
    using buffer_t = std::vector<ser_type>;
    
  public:

    NDVecArrayStreamer(std::size_t initial_buffer_size = 10000);

    void serialize(std::fstream& stream, const type& val, const shape_t& chunk_size, const StreamerMode& mode);
    void deserialize(std::fstream& stream, type& val);

    // appends slice to the end of the corresponding axis
    // only works if chunk_size is also a slice and the passed `chunk` has the same size
    void append_slice(std::fstream& stream, const type& chunk, const StreamerMode& mode);

  private:
    
    void serialize_all_chunks_dense(std::fstream& stream, const type& val, const shape_t& chunk_size);
    void serialize_all_chunks_zero_suppressed(std::fstream& stream, const type& val, const shape_t& chunk_size);

  private:

    void write_buffer(std::size_t num_elems, std::fstream& stream);
    void read_into_buffer(std::size_t num_elems, std::fstream& stream);
    
  private:

    std::shared_ptr<buffer_t> m_ser_buffer;

  };
}

#include "NDVecArrayStreamer.hxx"
