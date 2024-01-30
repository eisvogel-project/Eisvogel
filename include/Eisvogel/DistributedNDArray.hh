#ifndef __DISTRIBUTED_NDARRAY_HH
#define __DISTRIBUTED_NDARRAY_HH

#include <iostream>

#include <string>
#include <memory>
#include <map>
#include <queue>
#include <filesystem>
#include <fstream>

#include "Eisvogel/NDArray.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "Eisvogel/Serialization.hh"

struct ChunkMetadata {

  ChunkMetadata(const std::string filename, const IndexVector& start_ind, const IndexVector& stop_ind) :
    filename(filename), start_ind(start_ind), stop_ind(stop_ind) { } 
  
  std::string filename;
  IndexVector start_ind;
  IndexVector stop_ind;  
};

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

// ----

template <class T, std::size_t dims>
class DistributedNDArray : public NDArray<T, dims> {

public:
  
  DistributedNDArray(std::string dirpath, std::size_t max_cache_size);
  ~DistributedNDArray();

  using chunk_t = DenseNDArray<T, dims>;

  // For assembling a distributed array
  void RegisterChunk(const chunk_t& chunk, const IndexVector start_ind);
  void Flush();
  
  // For accessing a distributed array
  T operator()(IndexVector& inds);

private:

  bool chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds);
  std::size_t getChunkIndex(const IndexVector& inds);
  chunk_t& retrieveChunk(std::size_t chunk_ind);
  void calculateShape();
  bool isGloballyContiguous(IndexVector& global_start_inds, IndexVector& global_stop_inds);

  IndexVector& getGlobalStartInd();
  IndexVector& getGlobalStopInd();
  std::size_t getVolume(IndexVector& start_inds, IndexVector& stop_inds);
  
private:

  const std::string m_dirpath;
  const std::string m_indexpath;
  const std::size_t m_max_cache_size;

  // Keeps track of the chunks this DistributedNDArray is composed of
  using index_t = std::vector<ChunkMetadata>;
  index_t m_chunk_index;

  // Data strutures for caching of frequently-accessed elements of the array
  std::map<std::size_t, chunk_t> m_chunk_cache; // key is index of chunk in m_chunk_index
  std::queue<std::size_t> m_cache_queue; // to keep track of the age of cached chunks
};

#include "DistributedNDArray.hxx"

#endif
