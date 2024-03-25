#pragma once

#include "NDVecArray.hh"

namespace stor {

  enum class StreamerMode : uint32_t {
    automatic = 0,
    dense = 1,
    zero_suppressed = 2
  };
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  class NDVecArrayStreamer {

    using type = NDVecArray<T, dims, vec_dims>;
    using data_t = typename type::data_t;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;   

  public:

    NDVecArrayStreamer(const shape_t& chunk_size);

    static void serialize(std::fstream& stream, const type& val, const StreamerMode& mode);
    static void deserialize(std::fstream& stream, type& val);

    // only works, if chunk_size is also a slice and the passed `chunk` has the same size
    static void append_slice(std::fstream& stream, const type& chunk);

    // to be made private and called from functions above
    
    static void serialize_dense(std::fstream& stream, const type& val);
    static void deserialize_dense(std::fstream& stream, type& val);
    static type deserialize_dense(std::fstream& stream);

    // chunk_size for (de)serialization
    static void serialize_suppress_zero(std::fstream& stream, const type& val);
    static void deserialize_suppress_zero(std::fstream& stream, type& val);
    static type deserialize_suppress_zero(std::fstream& stream);

    // on-disk operations


  private:

    

  };
}

#include "NDVecArrayStreamer.hxx"
