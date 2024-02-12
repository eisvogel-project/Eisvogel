#include <uuid/uuid.h>
#include <type_traits>

namespace stor {
  template <>
  struct Traits<ChunkMetadata> {
    using type = ChunkMetadata;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<uint32_t>::serialize(stream, static_cast<uint32_t>(val.chunk_type));
      Traits<std::string>::serialize(stream, val.filename);
      Traits<IndexVector>::serialize(stream, val.start_ind);
      Traits<IndexVector>::serialize(stream, val.stop_ind);
    }

    static type deserialize(std::iostream& stream) {
      ChunkType chunk_type{Traits<uint32_t>::deserialize(stream)};
      std::string filename = Traits<std::string>::deserialize(stream);
      IndexVector start_ind = Traits<IndexVector>::deserialize(stream);
      IndexVector stop_ind = Traits<IndexVector>::deserialize(stream);
      return ChunkMetadata(filename, start_ind, stop_ind, chunk_type);
    }    
  };
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::DistributedNDArray(std::string dirpath, std::size_t max_cache_size, SerializerT& ser) :
  NDArray<T, dims>(), m_dirpath(dirpath), m_indexpath(dirpath + "/index.bin"), m_max_cache_size(max_cache_size),
  m_global_start_ind(dims, 0), m_ser(ser) {

  // Create directory if it does not already exist
  if(!std::filesystem::exists(m_dirpath)) {
    std::filesystem::create_directory(m_dirpath);
  }

  if(std::filesystem::exists(m_indexpath)) {
    // Load index if there is one
    std::fstream ifs; 
    ifs.open(m_indexpath, std::ios::in | std::ios::binary);
    m_chunk_index = m_ser.template deserialize<index_t>(ifs);
    ifs.close();    
  }
  else {
    // Attempt to rebuild index
    rebuildIndex();
  }

  // Calculate global shape of this array
  calculateShape();
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::~DistributedNDArray() { }

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
template <class ChunkT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::RegisterChunk(const ChunkT& chunk, const IndexVector start_ind, bool require_nonoverlapping) {

  IndexVector stop_ind = start_ind + chunk.shape();
  
  if(require_nonoverlapping) {
    // make sure this chunk does not overlap with any that we already have and crash if it does
    for(auto& chunk_meta : m_chunk_index) {
      if(chunkContainsInds(chunk_meta, start_ind) || chunkContainsInds(chunk_meta, stop_ind)) {
	throw std::runtime_error("Trying to add a chunk that overlaps with an already existing one!");
      }
    }
  }

  WriteChunk(chunk, start_ind, true);
  
  // update global shape of this array (if possible)
  calculateShape();
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
template <class ChunkT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::WriteChunk(const ChunkT& chunk, const IndexVector start_ind, bool add_to_index) {

  IndexVector stop_ind = start_ind + chunk.shape();
  
  // get chunk filename that should not clash with anything
  uuid_t uuid_binary;
  uuid_generate_random(uuid_binary);
  char uuid_string[36];
  uuid_unparse(uuid_binary, uuid_string);
  std::string chunk_filename = std::string(uuid_string) + ".bin";

  // build chunk metadata and add to index
  ChunkType chunk_type = ChunkType::dense;
  if constexpr(std::is_same_v<ChunkT, sparse_t>) {
    chunk_type = ChunkType::sparse;
  }
  ChunkMetadata meta(chunk_filename, start_ind, stop_ind, chunk_type);

  if(add_to_index) {
    m_chunk_index.push_back(meta);
  }

  // write data to disk
  std::string chunk_path = m_dirpath + "/" + chunk_filename;
  std::fstream ofs;
  ofs.open(chunk_path, std::ios::out | std::ios::binary);  
  m_ser.template serialize<ChunkMetadata>(ofs, meta);   
  m_ser.template serialize<ChunkT>(ofs, chunk);
  ofs.close();  
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::MakeIndexPersistent() {
  // Update index on disk
  std::fstream ofs;
  ofs.open(m_indexpath, std::ios::out | std::ios::binary);  
  m_ser.template serialize<index_t>(ofs, m_chunk_index);
  ofs.close();
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::rebuildIndex() {

  // clear the cached index
  m_chunk_index.clear();
  
  // if there is a persistent index file, make sure to delete it to prevent inconsistencies
  if(std::filesystem::exists(m_indexpath)) {
    std::filesystem::remove(m_indexpath);
  }

  // With the index gone, also the cache is now out of scope
  m_chunk_cache.clear();
  m_cache_queue = {};

  for(auto const& dir_entry : std::filesystem::directory_iterator(m_dirpath)) {
    std::fstream ifs;
    ifs.open(dir_entry.path(), std::ios::in | std::ios::binary);
    ChunkMetadata meta = m_ser.template deserialize<ChunkMetadata>(ifs);
    ifs.close();
    
    m_chunk_index.push_back(meta);
  }
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::RebuildChunks(const IndexVector& requested_chunk_size) {

  // How to make sure own index doesn't get messed up when new chunks are generated?

  // 0) factor out most things of currently existing RegisterChunk into new private method:
  //      -> WriteChunk(chunk, start_ind, add_to_index)
  //   (keep the overlap checking etc. in the original RegisterChunk)
  
  // Then, this function will look like
  // 0) rebuild index (to make sure everything is taken into account)
  // 1) then loop over chunks like in `WeightingFieldUtils`
  //     -> get chunk that has the starting indices
  //     -> if it conforms to the requested size, don't touch it and continue (will guarantee this function is idempotent)
  //     -> if it doesn't, build a chunk of correct size (using NDArray::range)
  // 2) write the individual chunks, but don't touch the index (using `WriteChunk`)
  // 3) remove all chunks in the index (-> this will remove all the old ones)

  std::cout << "in RebuildChunks" << std::endl;
  
  if(requested_chunk_size.size() != dims) {
    throw std::runtime_error("Error: requested chunk size has wrong dimensionality!");
  }

  std::cout << "rebuilding index" << std::endl;
  
  // Make sure we start from a clean index
  rebuildIndex(); 

  std::cout << "rebuilt index" << std::endl;
  
  if(!isGloballyContiguous(getGlobalStartInd(), getGlobalStopInd())) {
    throw std::runtime_error("Error: refusing to rebuild chunks for a non-contiguous array!");
  }

  calculateShape();

  IndexVector global_shape(this -> m_shape);

  std::cout << "global shape" << std::endl;
  global_shape.print();
  
  IndexVector number_required_chunks = (global_shape + requested_chunk_size - 1) / requested_chunk_size;

  std::cout << "will have " << std::endl;
  number_required_chunks.print();
  std::cout << " chunks after rebuilding" << std::endl;
  
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
T DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::operator()(IndexVector& inds) {

  // check to which chunk this index belongs
  std::size_t chunk_ind = getChunkIndex(inds);

  // retrieve chunk from cache or load from file
  dense_t& found_chunk = retrieveChunk(chunk_ind);
  
  // index and return element
  IndexVector inds_within_chunk = inds - m_chunk_index[chunk_ind].start_ind;  
  return found_chunk(inds_within_chunk);
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
bool DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds) {
  return isInIndexRange(inds, chunk_meta.start_ind, chunk_meta.stop_ind);
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
std::size_t DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getChunkIndex(const IndexVector& inds) {
  std::size_t chunk_ind = 0;
  for(chunk_ind = 0; chunk_ind < m_chunk_index.size(); chunk_ind++) {
    if(chunkContainsInds(m_chunk_index[chunk_ind], inds)) {
      return chunk_ind;
    }
  }
  
  throw std::runtime_error("No chunk provides these indices!");
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::dense_t& DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::retrieveChunk(std::size_t chunk_ind) {
  
  ChunkMetadata& chunk_meta = m_chunk_index[chunk_ind];

  if(!m_chunk_cache.contains(chunk_ind)) {
    
    // do we need to free some space in the cache?
    if(m_chunk_cache.size() > m_max_cache_size) {
      std::size_t dense_to_free = m_cache_queue.front();
      std::cout << "Removing chunk " + std::to_string(dense_to_free) + " from cache" << std::endl;
      m_chunk_cache.erase(dense_to_free);
      m_cache_queue.pop();
    }    
    
    std::fstream ifs;
    std::string chunk_path = m_dirpath + "/" + chunk_meta.filename;
    std::cout << "Loading chunk from " + chunk_path + " ... ";
    ifs.open(chunk_path, std::ios::in | std::ios::binary);
    ChunkMetadata meta = m_ser.template deserialize<ChunkMetadata>(ifs);

    if(meta.chunk_type == ChunkType::dense) {
      m_chunk_cache.insert({chunk_ind, m_ser.template deserialize<dense_t>(ifs)});
    }
    else if(meta.chunk_type == ChunkType::sparse) {
      m_chunk_cache.insert({chunk_ind, dense_t::From(m_ser.template deserialize<sparse_t>(ifs))});
    }
    else {
      throw std::runtime_error("Error: unknown chunk type encountered!");
    }
    
    m_cache_queue.push(chunk_ind);
    ifs.close();
    std::cout << "done!" << std::endl;
  }
  
  return m_chunk_cache.find(chunk_ind) -> second;
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::calculateShape() {

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

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
bool DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::isGloballyContiguous(IndexVector& global_start_inds, IndexVector& global_stop_inds) {
  
  std::size_t global_volume = getVolume(global_start_inds, global_stop_inds);
  
  std::size_t total_chunk_volume = 0;
  for(auto& chunk : m_chunk_index) {
    total_chunk_volume += getVolume(chunk.start_ind, chunk.stop_ind);
  }

  return global_volume == total_chunk_volume;
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
std::size_t DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getVolume(IndexVector& start_inds, IndexVector& stop_inds) {
  std::size_t volume = 1;
  for(std::size_t cur_dim = 0; cur_dim < dims; cur_dim++) {
    volume *= (stop_inds(cur_dim) - start_inds(cur_dim));
  }
  return volume;
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
IndexVector& DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getGlobalStartInd() {
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

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
IndexVector& DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getGlobalStopInd() {
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
