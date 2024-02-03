#include <uuid/uuid.h>

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
DistributedNDArray<T, dims>::DistributedNDArray(std::string dirpath, std::size_t max_cache_size) :
  NDArray<T, dims>(), m_dirpath(dirpath), m_indexpath(dirpath + "/index.bin"), m_max_cache_size(max_cache_size),
  m_global_start_ind(dims, 0) {

  // Create directory if it does not already exist
  if(!std::filesystem::exists(m_dirpath)) {
    std::filesystem::create_directory(m_dirpath);
  }

  if(std::filesystem::exists(m_indexpath)) {
    // Load index if there is one
    std::fstream ifs; 
    ifs.open(m_indexpath, std::ios::in | std::ios::binary);
    stor::Serializer iser(ifs);
    m_chunk_index = iser.deserialize<index_t>();
    ifs.close();    
  }
  else {
    // Attempt to rebuild index
    rebuildIndex();
  }

  // Calculate global shape of this array
  calculateShape();
}

template <class T, std::size_t dims>
DistributedNDArray<T, dims>::~DistributedNDArray() { }

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::RegisterChunk(const DenseNDArray<T, dims>& chunk, const IndexVector start_ind, bool require_nonoverlapping) {

  IndexVector stop_ind = start_ind + chunk.shape();
  
  if(require_nonoverlapping) {
    // make sure this chunk does not overlap with any that we already have and crash if it does
    for(auto& chunk_meta : m_chunk_index) {
      if(chunkContainsInds(chunk_meta, start_ind) || chunkContainsInds(chunk_meta, stop_ind)) {
	throw std::runtime_error("Trying to add a chunk that overlaps with an already existing one!");
      }
    }
  }

  // get chunk filename that should not clash with anything
  uuid_t uuid_binary;
  uuid_generate_random(uuid_binary);
  char uuid_string[36];
  uuid_unparse(uuid_binary, uuid_string);
  std::string chunk_filename = std::string(uuid_string) + ".bin";

  // build chunk metadata and add to index
  ChunkMetadata meta(chunk_filename, start_ind, stop_ind);
  m_chunk_index.push_back(meta);
  
  // write chunk data to disk --> this is where fancy sparsification and compression would happen
  
  std::string chunk_path = m_dirpath + "/" + chunk_filename;
  std::fstream ofs;
  ofs.open(chunk_path, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);
  oser.serialize<ChunkMetadata>(meta);
  oser.serialize<chunk_t>(chunk);
  ofs.close();

  // update global shape of this array (if possible)
  calculateShape();
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::MakeIndexPersistent() {
  // Update index on disk
  std::fstream ofs;
  ofs.open(m_indexpath, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);
  oser.serialize<index_t>(m_chunk_index);
  ofs.close();
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::rebuildIndex() {
  m_chunk_index.clear();

  // With the index gone, also the cache is now out of scope
  m_chunk_cache.clear();
  m_cache_queue = {};

  for(auto const& dir_entry : std::filesystem::directory_iterator(m_dirpath)) {
    std::fstream ifs;
    ifs.open(dir_entry.path(), std::ios::in | std::ios::binary);
    stor::Serializer iser(ifs);
    ChunkMetadata meta = iser.deserialize<ChunkMetadata>();
    ifs.close();
    
    m_chunk_index.push_back(meta);
  }
}

template <class T, std::size_t dims>
T DistributedNDArray<T, dims>::operator()(IndexVector& inds) {

  // check to which chunk this index belongs
  std::size_t chunk_ind = getChunkIndex(inds);

  // retrieve chunk from cache or load from file
  chunk_t& found_chunk = retrieveChunk(chunk_ind);
  
  // index and return element
  IndexVector inds_within_chunk = inds - m_chunk_index[chunk_ind].start_ind;  
  return found_chunk(inds_within_chunk);
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
    std::cout << "Loading chunk from " + chunk_path + " ... ";
    ifs.open(chunk_path, std::ios::in | std::ios::binary);
    stor::Serializer iser(ifs);
    iser.deserialize<ChunkMetadata>(); // skip metadata
    m_chunk_cache.insert({chunk_ind, iser.deserialize<chunk_t>()});
    m_cache_queue.push(chunk_ind);
    ifs.close();
    std::cout << "done!" << std::endl;
  }
  
  return m_chunk_cache.find(chunk_ind) -> second;
}

template <class T, std::size_t dims>
void DistributedNDArray<T, dims>::calculateShape() {

  if(m_chunk_index.empty()) {
    // Nothing to do if everything is empty
    return;
  }
  
  m_global_start_ind = getGlobalStartInd();
  IndexVector& global_stop_ind = getGlobalStopInd();

  if(isGloballyContiguous(m_global_start_ind, global_stop_ind)) {
    // Chunks fill a contiguous array, makes sense to define a global shape
    for(std::size_t cur_dim = 0; cur_dim < dims; cur_dim++) {
      this -> m_shape[cur_dim] = global_stop_ind(cur_dim) - m_global_start_ind(cur_dim);
    }
  }
}

template <class T, std::size_t dims>
bool DistributedNDArray<T, dims>::isGloballyContiguous(IndexVector& global_start_inds, IndexVector& global_stop_inds) {
  
  std::size_t global_volume = getVolume(global_start_inds, global_stop_inds);
  
  std::size_t total_chunk_volume = 0;
  for(auto& chunk : m_chunk_index) {
    total_chunk_volume += getVolume(chunk.start_ind, chunk.stop_ind);
  }

  return global_volume == total_chunk_volume;
}

template <class T, std::size_t dims>
std::size_t DistributedNDArray<T, dims>::getVolume(IndexVector& start_inds, IndexVector& stop_inds) {
  std::size_t volume = 1;
  for(std::size_t cur_dim = 0; cur_dim < dims; cur_dim++) {
    volume *= (stop_inds(cur_dim) - start_inds(cur_dim));
  }
  return volume;
}

template <class T, std::size_t dims>
IndexVector& DistributedNDArray<T, dims>::getGlobalStartInd() {
  if(m_chunk_index.empty()) {
    throw std::runtime_error("Trying to compute start index of an empty array!");
  }
  IndexVector* global_start_ind = &m_chunk_index[0].start_ind;
  for(auto& chunk : m_chunk_index) {
    if(all(chunk.start_ind <= *global_start_ind)) {
      global_start_ind = &chunk.start_ind;
    }
  }
  return *global_start_ind;
}

template <class T, std::size_t dims>
IndexVector& DistributedNDArray<T, dims>::getGlobalStopInd() {
  if(m_chunk_index.empty()) {
    throw std::runtime_error("Trying to compute stop index of an empty array!");
  }
  IndexVector* global_stop_ind = &m_chunk_index[0].stop_ind;
  for(auto& chunk : m_chunk_index) {
    if(all(chunk.stop_ind >= *global_stop_ind)) {
      global_stop_ind = &chunk.stop_ind;
    }
  }
  return *global_stop_ind;
}
