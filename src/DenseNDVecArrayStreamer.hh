#pragma once

#include "DenseNDVecArray.hh"

namespace stor {

  // has some specialized functionality for efficient serialization / deserialization / on-disk manipulations
  // that the normal DefaultSerializer does not have

  // TODO: experimental for now, should put into different structs for dense / suppress_zero / etc.?
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  class DenseNDVecArrayStreamer {

    using type = DenseNDVecArray<T, dims, vec_dims>;
    using data_t = typename type::data_t;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;   

  public:

    DenseNDVecArrayStreamer(const shape_t& chunk_size);
    
    static void serialize_dense(std::fstream& stream, const type& val);
    static void deserialize_dense(std::fstream& stream, type& val);
    static type deserialize_dense(std::fstream& stream);

    // chunk_size for (de)serialization
    static void serialize_suppress_zero(std::fstream& stream, const type& val);
    static void deserialize_suppress_zero(std::fstream& stream, type& val);
    static type deserialize_suppress_zero(std::fstream& stream);

    // on-disk operations
    static void append_chunk(std::fstream& stream, const type& chunk);

  private:
   
    // (de)serialization buffers
  };
}

#include "DenseNDVecArrayStreamer.hxx"
