#ifndef __DISTRIBUTED_NDARRAY_HH
#define __DISTRIBUTED_NDARRAY_HH

#include <iostream>

#include <string>
#include <memory>
#include <map>
#include <queue>
#include <filesystem>
#include <fstream>

#include "DenseNDArray.hh"
#include "SparseNDArray.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "Serialization.hh"

// global TODO: remove SparseT and DenseT and instead have a generic type argument for the to-be-saved chunk
enum class ChunkType : uint32_t {
  dense = 0,
  sparse = 1
};

struct ChunkNotFoundError : public std::runtime_error {
  ChunkNotFoundError() : std::runtime_error("No chunk provides these indices!") { };
};

struct ChunkMetadata {

  ChunkMetadata(const std::string filename, const IndexVector& start_ind, const IndexVector& stop_ind, const ChunkType& chunk_type) :
    filename(filename), start_ind(start_ind), stop_ind(stop_ind), chunk_type(chunk_type) { } 

  bool operator==(const ChunkMetadata& rhs) {
    if(chunk_type != rhs.chunk_type) return false;
    if(filename != rhs.filename) return false;
    if(start_ind != rhs.start_ind) return false;
    if(stop_ind != rhs.stop_ind) return false;
    return true;
  }
  
  std::string filename;
  IndexVector start_ind;
  IndexVector stop_ind;
  ChunkType chunk_type;
};

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template <class, std::size_t> class SparseT, class SerializerT>
class DistributedNDArray : public NDArray<T, dims> {

public:

  using dense_t = DenseT<T, dims>;
  using sparse_t = SparseT<T, dims>;
  
  DistributedNDArray(std::string dirpath, std::size_t max_cache_size, SerializerT& ser);
  ~DistributedNDArray();

  // For assembling and indexing a distributed array
  template <class ChunkT>
  void RegisterChunk(const ChunkT& chunk, const IndexVector start_ind, bool require_nonoverlapping = false);    
  void MakeIndexPersistent();
  void rebuildIndex();
  
  // For accessing a single element
  T operator()(IndexVector& inds);

  // For accessing a range of elements
  dense_t range(const IndexVector& start_inds, const IndexVector& stop_inds);
  
  std::size_t startInd(std::size_t dim) const {
    return m_global_start_ind(dim);
  }

  // Attempts to redistribute into equally large chunks
  void RebuildChunks(const IndexVector& requested_chunk_size);

  // Simpler version that only merges chunks
  void MergeChunks(std::size_t dim, std::size_t max_dimsize);

  void printChunks();
  
private:

  template <class ChunkT>
  void WriteChunk(const ChunkT& chunk, const IndexVector start_ind, bool add_to_index = true);
  
  bool chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds);
  std::size_t getChunkIndex(const IndexVector& inds);
  std::size_t getNeighbouringChunkIndex(std::size_t chunk_index, std::size_t dim);
  dense_t& retrieveChunk(std::size_t chunk_ind);
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
  // TODO: maybe later when we need fancier things (e.g. predictive loading of additional neighbouring chunks),
  // can think about turning this into a class
  index_t m_chunk_index;

  // The index may not start at {0, 0, 0}
  IndexVector m_global_start_ind;
  
  // Data strutures for caching of frequently-accessed elements of the array
  std::map<std::size_t, dense_t> m_chunk_cache; // key is index of chunk in m_chunk_index
  std::queue<std::size_t> m_cache_queue; // to keep track of the age of cached chunks

  SerializerT& m_ser;
};

#include "DistributedNDArray.hxx"

template <class T, std::size_t dims, class SerializerT = stor::DefaultSerializer>
using DistributedScalarNDArray = DistributedNDArray<T, dims, DenseNDArray, SparseNDArray, SerializerT>;

#endif
