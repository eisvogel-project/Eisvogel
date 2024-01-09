#ifndef __DISTRIBUTED_NDARRAY_HH
#define __DISTRIBUTED_NDARRAY_HH

#include <iostream>

#include <string>
#include <memory>
#include <map>

#include "Eisvogel/NDArray.hh"
#include "Eisvogel/Serialization.hh"

struct ChunkMetadata {

  ChunkMetadata(const std::string filename, const IndexVector& start_ind, const IndexVector& stop_ind) :
    filename(filename), start_ind(start_ind), stop_ind(stop_ind) { } 
  
  std::string filename;
  IndexVector start_ind;
  IndexVector stop_ind;  
};

// ----

namespace stor {
  template <>
  struct Traits<ChunkMetadata> {
    using type = ChunkMetadata;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<std::string>::serialize(stream, val.filename);
      Traits<IndexVector>::serialize(stream, val.start_ind);
      Traits<IndexVector>::serialize(stream, val.stop_ind);
    }

    static type deserialize(std::iostream& stream) {
      std::string filename = Traits<std::string>::deserialize(stream);
      IndexVector start_ind = Traits<IndexVector>::deserialize(stream);
      IndexVector stop_ind = Traits<IndexVector>::deserialize(stream);
      return ChunkMetadata(filename, start_ind, stop_ind);
    }    
  };
}

template <class T, std::size_t dims>
class DistributedNDArray {

public:
  
  DistributedNDArray(std::string dirpath, std::size_t max_cache_size);

  using chunk_t = NDArray<T, dims>;

  // For assembling a distributed array
  void RegisterChunk(const chunk_t& chunk, const IndexVector start_ind);
  
  // For accessing a distributed array
  T& operator()(DenseVector<T>& inds);

private:

  std::string m_dirpath;
  std::size_t m_max_cache_size;

  std::vector<ChunkMetadata> m_chunk_index;

  using chunk_cache_t = DenseNDArray<T, dims>;
  std::map<std::size_t, chunk_cache_t> m_chunk_cache; // key is index of chunk in m_chunk_index    
};

// ---

template <class T, std::size_t dims>
DistributedNDArray<T, dims>::DistributedNDArray(std::string dirpath, std::size_t max_cache_size) :
  m_dirpath(dirpath), m_max_cache_size(max_cache_size) {

  // create directory if it does not exist
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::RegisterChunk(const NDArray<T, dims>& chunk, const IndexVector start_ind) {

  // make sure this chunk does not overlap with any that we already have

  // add chunk metadata to index

  // write chunk data
  
}

template <class T, std::size_t dims>
T& DistributedNDArray<T, dims>::operator()(DenseVector<T>& inds) {

  // check to which chunk this index belongs

  // check if chunk is in cache

  // retrieve chunk from cache or load from file

  // index and return element
  
}

#endif
