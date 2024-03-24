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
    
    static void serialize_dense(std::fstream& stream, const type& val);  
    static type deserialize_dense(std::fstream& stream);
    
    static void serialize_suppress_zero(std::fstream& stream, const type& val, const Vector<std::size_t, dims>& chunk_size);
    static type deserialize_suppress_zero(std::fstream& stream, const Vector<std::size_t, dims>& chunk_size);

  private:
   
    // buffers
  };
}

#include "DenseNDVecArrayStreamer.hxx"
