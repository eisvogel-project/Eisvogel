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
class DistributedNDArray {

public:
  
  DistributedNDArray(std::string dirpath, std::size_t max_cache_size);
  ~DistributedNDArray();

  using chunk_t = DenseNDArray<T, dims>;
  using shape_t = std::array<std::size_t, dims>;

  // For assembling a distributed array
  void RegisterChunk(const chunk_t& chunk, const IndexVector start_ind);
  void Flush();
  
  // For accessing a distributed array
  T& operator()(IndexVector& inds);
  const shape_t& shape();

private:

  bool chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds);
  std::size_t getChunkIndex(const IndexVector& inds);
  chunk_t& retrieveChunk(std::size_t chunk_ind);
  void calculateShape();
  bool isContiguous();

protected:
  const shape_t m_shape{};
  
private:

  const std::string m_dirpath;
  const std::string m_indexpath;
  const std::size_t m_max_cache_size;

  using index_t = std::vector<ChunkMetadata>;
  index_t m_chunk_index;

  std::map<std::size_t, chunk_t> m_chunk_cache; // key is index of chunk in m_chunk_index
  std::queue<std::size_t> m_cache_queue; // to keep track of the age of cached chunks
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

  // make sure this chunk does not overlap with any that we already have and crash if it does
  IndexVector stop_ind = start_ind + chunk.shape();
  for(auto& chunk_meta : m_chunk_index) {
    if(chunkContainsInds(chunk_meta, start_ind) || chunkContainsInds(chunk_meta, stop_ind)) {
      throw std::runtime_error("Trying to add a chunk that overlaps with an already existing one!");
    }
  }

  // build chunk metadata and add to index
  std::string chunk_filename = "chunk_" + std::to_string(m_chunk_index.size()) + ".bin";
  ChunkMetadata meta(chunk_filename, start_ind, stop_ind);
  m_chunk_index.push_back(meta);
  
  // write chunk data to disk --> this is where fancy sparsification and compression would happen
  std::string chunk_path = m_dirpath + "/" + chunk_filename;
  std::fstream ofs;
  ofs.open(chunk_path, std::ios::out | std::ios::binary);  
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
T& DistributedNDArray<T, dims>::operator()(IndexVector& inds) {

  // check to which chunk this index belongs
  std::size_t chunk_ind = getChunkIndex(inds);

  // retrieve chunk from cache or load from file
  chunk_t& found_chunk = retrieveChunk(chunk_ind);

  std::cout << "found element in chunk " << std::to_string(chunk_ind) << std::endl;
  
  // index and return element
  return found_chunk(inds);
}

template <class T, std::size_t dims>
bool DistributedNDArray<T, dims>::chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds) {
  return isInIndexRange(inds, chunk_meta.start_ind, chunk_meta.stop_ind);
}

template <class T, std::size_t dims>
std::size_t DistributedNDArray<T, dims>::getChunkIndex(const IndexVector& inds) {
  std::size_t chunk_ind = 0;
  for(chunk_ind = 0; chunk_ind < m_chunk_index.size(); chunk_ind++) {
    if(chunkContainsInds(m_chunk_index[chunk_ind], inds)) {
      return chunk_ind;
    }
  }
  
  throw std::runtime_error("No chunk provides these indices!");
}

template <class T, std::size_t dims>
DistributedNDArray<T, dims>::chunk_t& DistributedNDArray<T, dims>::retrieveChunk(std::size_t chunk_ind) {

  ChunkMetadata& chunk_meta = m_chunk_index[chunk_ind];

  if(!m_chunk_cache.contains(chunk_ind)) {
    
    // do we need to free some space in the cache?
    if(m_chunk_cache.size() > m_max_cache_size) {
      std::size_t chunk_to_free = m_cache_queue.front();
      std::cout << "Removing chunk " + std::to_string(chunk_to_free) + " from cache" << std::endl;
      m_chunk_cache.erase(chunk_to_free);
      m_cache_queue.pop();
    }    
    
    // load chunk from file --> this is where fancy re-densification and decompression would happen
    std::fstream ifs;
    std::string chunk_path = m_dirpath + "/" + chunk_meta.filename;
    std::cout << "Loading chunk from " + chunk_path << std::endl;
    ifs.open(chunk_path, std::ios::in | std::ios::binary);
    stor::Serializer iser(ifs);
    m_chunk_cache.insert({chunk_ind, iser.deserialize<chunk_t>()});
    m_cache_queue.push(chunk_ind);
    ifs.close();    
  }
  
  return m_chunk_cache.find(chunk_ind) -> second;
}

template <class T, std::size_t dims>
const DistributedNDArray<T, dims>::shape_t& DistributedNDArray<T, dims>::shape() {
  return m_shape;
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::calculateShape() {

  // go through chunks and determine global start and end inds
  
}

template <class T, std::size_t dims>
bool DistributedNDArray<T, dims>::isContiguous() {

  // determine global (distributed-array wide) start and stop inds
  // determine total number of elements
  // check if sum of elements of chunks give teh same value  
}

#endif
