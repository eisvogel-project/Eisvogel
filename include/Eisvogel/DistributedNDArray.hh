#ifndef __DISTRIBUTED_NDARRAY_HH
#define __DISTRIBUTED_NDARRAY_HH

#include <iostream>

#include <string>
#include <memory>
#include <map>
#include <filesystem>
#include <fstream>

#include "Eisvogel/NDArray.hh"
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
class DistributedNDArray {

public:
  
  DistributedNDArray(std::string dirpath, std::size_t max_cache_size);
  ~DistributedNDArray();

  using chunk_t = DenseNDArray<T, dims>;

  // For assembling a distributed array
  void RegisterChunk(const chunk_t& chunk, const IndexVector start_ind);
  void Flush();
  
  // For accessing a distributed array
  T& operator()(DenseVector<T>& inds);

private:

  const std::string m_dirpath;
  const std::string m_indexpath;
  const std::size_t m_max_cache_size;

  using index_t = std::vector<ChunkMetadata>;
  index_t m_chunk_index;

  std::map<std::size_t, chunk_t> m_chunk_cache; // key is index of chunk in m_chunk_index    
};

// ---

template <class T, std::size_t dims>
DistributedNDArray<T, dims>::DistributedNDArray(std::string dirpath, std::size_t max_cache_size) :
  m_dirpath(dirpath), m_indexpath(dirpath + "/index.bin"), m_max_cache_size(max_cache_size) {

  // Create directory if it does not already exist
  if(!std::filesystem::exists(m_dirpath)) {
    std::filesystem::create_directory(m_dirpath);
  }

  // Load index if there is one
  if(std::filesystem::exists(m_indexpath)) {
    // Load index
    std::fstream ifs; 
    ifs.open(m_indexpath, std::ios::in | std::ios::binary);
    stor::Serializer iser(ifs);
    m_chunk_index = iser.deserialize<index_t>();
    ifs.close();    
  }
}

template <class T, std::size_t dims>
DistributedNDArray<T, dims>::~DistributedNDArray() {
  Flush();
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::RegisterChunk(const DenseNDArray<T, dims>& chunk, const IndexVector start_ind) {

  // TODO: make sure this chunk does not overlap with any that we already have and crash if it does

  // build chunk metadata and add to index
  std::string chunk_filename = m_dirpath + "/chunk_" + std::to_string(m_chunk_index.size()) + ".bin";
  IndexVector stop_ind = start_ind + chunk.shape();
  ChunkMetadata meta(chunk_filename, start_ind, stop_ind);
  m_chunk_index.push_back(meta);

  // write chunk data to disk
  std::fstream ofs;
  ofs.open(chunk_filename, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);
  oser.serialize<chunk_t>(chunk);
  ofs.close();  
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::Flush() {
  // Update index on disk
  std::fstream ofs;
  ofs.open(m_indexpath, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);
  oser.serialize<index_t>(m_chunk_index);
  ofs.close();
}

template <class T, std::size_t dims>
T& DistributedNDArray<T, dims>::operator()(DenseVector<T>& inds) {

  // check to which chunk this index belongs

  // check if chunk is in cache

  // retrieve chunk from cache or load from file

  // index and return element
  
}

#endif
