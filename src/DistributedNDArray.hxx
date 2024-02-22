#include <uuid/uuid.h>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "NDArrayOperations.hh"

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
  m_chunk_last_accessed(0), m_global_start_ind(dims, 0), m_ser(ser) {

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

  // get chunk filename that should not clash with anything
  uuid_t uuid_binary;
  uuid_generate_random(uuid_binary);
  char uuid_string[36];
  uuid_unparse(uuid_binary, uuid_string);
  std::string chunk_filename = std::string(uuid_string) + ".bin";
  
  IndexVector stop_ind = start_ind + chunk.shape();

  // prepare chunk metadata
  ChunkMetadata meta(chunk_filename, start_ind, stop_ind, ChunkType::dense);

  // prepare output file
  std::string chunk_path = m_dirpath + "/" + chunk_filename;
  std::fstream ofs;
  ofs.open(chunk_path, std::ios::out | std::ios::binary);  
  
  // Check number of nonzero elements to see if we should store this as dense or sparse array
  // TODO: later, check possibility of embellishing the number count onto the DenseArray to avoid recomputations
  std::size_t num_nonzero_elems = NDArrayOps::number_nonzero_elems(chunk);
  std::size_t num_elems = chunk.volume();

  // pays off to store as sparse chunk
  std::size_t sparse_vs_dense_expense_ratio = 3; // sparse storage is approximately 3x as expensive as dense storage per nonzero element
  if(sparse_vs_dense_expense_ratio * num_nonzero_elems < num_elems) {
    
    std::cout << "going to sparsify" << std::endl;
    
    auto to_keep = [](scalar_t value) -> bool {
      return value != 0.0;
    };        
    sparse_t sparse_chunk = sparse_t::From(chunk, to_keep, 0.0);    
    meta.chunk_type = ChunkType::sparse;

    std::cout << "after sparsification, " << sparse_chunk.NumEntries() << " entries remain" << std::endl;
    
    m_ser.template serialize<ChunkMetadata>(ofs, meta);   
    m_ser.template serialize<sparse_t>(ofs, sparse_chunk);
  }
  else {  
    std::cout << "store as dense" << std::endl;
    m_ser.template serialize<ChunkMetadata>(ofs, meta);   
    m_ser.template serialize<dense_t>(ofs, chunk);
  }
  
  if(add_to_index) {
    m_chunk_index.push_back(meta);
  }
  
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
  // 3) remove all chunks in the index (-> this will remove all the old ones), BUT not those that haven't been touched because they are already conforming to the correct size

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

  index_t chunks_to_keep;
  
  IndexVector buf_inds_start(dims, 0);
  IndexVector buf_inds_end = number_required_chunks;
  for(IndexCounter cnt(buf_inds_start, buf_inds_end); cnt.running(); ++cnt) {
    
    IndexVector chunk_ind = cnt.index();
    IndexVector chunk_inds_start = chunk_ind * requested_chunk_size;
    IndexVector chunk_inds_end = NDArrayOps::min(chunk_inds_start + requested_chunk_size, global_shape);
    
    IndexVector actual_chunk_shape = chunk_inds_end - chunk_inds_start;

    // check if the current chunk already conforms to the requirements
    std::size_t chunk_index = getChunkIndex(chunk_inds_start);
    dense_t current_chunk = retrieveChunk(chunk_index);

    if(actual_chunk_shape == current_chunk.shape()) {
      // std::cout << "chunk already has the correct size, keep it" << std::endl;
      chunks_to_keep.push_back(m_chunk_index[chunk_index]);
      continue;
    }
    
    // std::cout << "now working on rebuild chunk with inds" << std::endl;
    // std::cout << "chunk_inds_start = " << std::endl;
    chunk_inds_start.print();
    // std::cout << "chunk_inds_end = " << std::endl;
    chunk_inds_end.print();
    
    dense_t chunk = range(chunk_inds_start, chunk_inds_end);
    WriteChunk(chunk, chunk_inds_start, false);
  }

  std::cout << "now cleaning up old chunks" << std::endl;
  for(ChunkMetadata& to_keep : chunks_to_keep) {
    std::erase(m_chunk_index, to_keep);
  }
    
  for(ChunkMetadata& cur_meta : m_chunk_index) {
    std::string chunk_path = m_dirpath + "/" + cur_meta.filename;
    std::cout << "now removing " << chunk_path << std::endl;
    std::filesystem::remove(chunk_path);
  }

  // re-index again
  rebuildIndex();
  calculateShape();
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::MergeChunks(std::size_t dim_to_merge, std::size_t max_dimsize) {

  // 0) Rebuild index (to make sure we have the full picture)

  // 2) Then, start with chunk that has the global start inds

  // --> operate entirely on the index and prepare a full list of chunk mergings that need to happen (without actually doing anything)
  //     --> remove any elements from the index that are already taken care of
  // --> check if there are any problems; if not, go ahead and implement the chunk mergings
  // --> go to the chunk that neighbours the current main one and continue
  
  // 3) again, have the to_keep mechanism
  
  std::cout << "in MergeNeighbouringChunks" << std::endl;

  rebuildIndex(); 

  std::cout << "rebuilt index" << std::endl;
  
  if(!isGloballyContiguous(getGlobalStartInd(), getGlobalStopInd())) {
    throw std::runtime_error("Error: refusing to merge chunks for a non-contiguous array!");
  }

  calculateShape();

  // put chunks in order along the merging axis
  std::vector<std::size_t> chunk_order(m_chunk_index.size());
  std::iota(chunk_order.begin(), chunk_order.end(), 0);
  std::stable_sort(chunk_order.begin(), chunk_order.end(),
		   [&m_chunk_index = m_chunk_index, &dim_to_merge = dim_to_merge](std::size_t ind_1, std::size_t ind_2) {
		     return m_chunk_index[ind_1].start_ind(dim_to_merge) > m_chunk_index[ind_2].start_ind(dim_to_merge);
		   });

  std::cout << "chunks after sorting:" << std::endl;
  for(std::size_t chunk_ind : chunk_order) {
    m_chunk_index[chunk_ind].start_ind.print();
  }        

  index_t chunks_to_keep;
  
  // Now operate on the ordered list of chunks until everything is done
  while(chunk_order.size() > 0) {
    std::size_t cur_chunk_index = chunk_order.back();
    ChunkMetadata& cur_chunk_meta = m_chunk_index[cur_chunk_index];  
    IndexVector chunk_size = cur_chunk_meta.stop_ind - cur_chunk_meta.start_ind;

    // after we're done, this chunk will conform to all requirements, so remove it from future consideration
    std::erase(chunk_order, cur_chunk_index);
    
    // this chunk is already large enough, we're done with it
    if(chunk_size(dim_to_merge) >= max_dimsize) {
      chunks_to_keep.push_back(cur_chunk_meta);
      continue;
    }

    // chunk is not yet large enough, merge with additional ones
    dense_t output_chunk = retrieveChunk(cur_chunk_index);
    std::size_t neighbour_chunk_index = cur_chunk_index;
    while(true) {
      try {
	neighbour_chunk_index = getNeighbouringChunkIndex(neighbour_chunk_index, dim_to_merge);	
	std::erase(chunk_order, neighbour_chunk_index);

	std::cout << "merging chunk " << cur_chunk_index << " with chunk " << neighbour_chunk_index << std::endl;
	  
	output_chunk = NDArrayOps::concatenate(output_chunk, retrieveChunk(neighbour_chunk_index), dim_to_merge);

	// chunk is now large enough
	if(output_chunk.shape(dim_to_merge) >= max_dimsize) {
	  break;
	}	
      } catch(const ChunkNotFoundError& e) {
	// there are no more neighbours in this direction; stop here
	break;
      }
    }

    // write it to disk
    WriteChunk(output_chunk, cur_chunk_meta.start_ind, false);
  }

  // remove all the chunks in the old index that should not be kept around
  std::cout << "now cleaning up old chunks" << std::endl;

  for(ChunkMetadata& to_keep : chunks_to_keep) {
    std::erase(m_chunk_index, to_keep);
  }
  
  for(ChunkMetadata& cur_meta : m_chunk_index) {
    std::string chunk_path = m_dirpath + "/" + cur_meta.filename;
    std::cout << "now removing " << chunk_path << std::endl;
    std::filesystem::remove(chunk_path);
  }

  // re-index again
  rebuildIndex();
  calculateShape();
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
void DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::printChunks() {

  for(ChunkMetadata& cur_meta : m_chunk_index) {
    std::cout << "-----------------------" << std::endl;
    std::cout << "chunk: " << cur_meta.filename << std::endl;
    std::cout << "start_ind: " << std::endl;
    cur_meta.start_ind.print();
    std::cout << "stop_ind: " << std::endl;
    cur_meta.stop_ind.print();
    std::cout << "chunk_type: " << static_cast<uint32_t>(cur_meta.chunk_type) << std::endl;
    std::cout << "-----------------------" << std::endl;
  }  
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
std::size_t DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getNeighbouringChunkIndex(std::size_t chunk_index, std::size_t dim) {

  ChunkMetadata& chunk_meta = m_chunk_index[chunk_index];  
  IndexVector chunk_size_along_dim(dims, 0);
  chunk_size_along_dim(dim) = chunk_meta.stop_ind(dim) - chunk_meta.start_ind(dim);
  IndexVector neighbour_chunk_start_ind = chunk_meta.start_ind + chunk_size_along_dim;
  
  return getChunkIndex(neighbour_chunk_start_ind);
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
DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::dense_t DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::range(const IndexVector& start_inds, const IndexVector& stop_inds) {

  // TODO: very simple and very slow---only temporary for now and to be made faster by intelligently splicing together individual chunks
  
  if((start_inds.size() != dims) || (stop_inds.size() != dims)) {
    throw std::runtime_error("Error: not a possible range!");
  }
  
  IndexVector range_shape = stop_inds - start_inds;
  
  // --- this is just a crutch for now until we have fixed-size vectors
  std::array<std::size_t, dims> range_shape_crutch;
  std::copy(std::begin(range_shape), std::end(range_shape), std::begin(range_shape_crutch));
  // ---------
  
  DenseNDArray<T, dims> retval(range_shape_crutch, T());    
  for(IndexCounter cnt(start_inds, stop_inds); cnt.running(); ++cnt) {
    IndexVector cur_ind = cnt.index();
    IndexVector range_ind = cur_ind - start_inds;
    retval(range_ind) = this -> operator()(cur_ind);
  }
  
  return retval;
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
bool DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::chunkContainsInds(const ChunkMetadata& chunk_meta, const IndexVector& inds) {
  return isInIndexRange(inds, chunk_meta.start_ind, chunk_meta.stop_ind);
}

template <class T, std::size_t dims, template<class, std::size_t> class DenseT, template<class, std::size_t> class SparseT, class SerializerT>
std::size_t DistributedNDArray<T, dims, DenseT, SparseT, SerializerT>::getChunkIndex(const IndexVector& inds) {

  if(m_chunk_index.size() == 0) {
    [[unlikely]];
    throw ChunkNotFoundError();
  }
  
  if(chunkContainsInds(m_chunk_index[m_chunk_last_accessed], inds)) {
    [[likely]];
    return m_chunk_last_accessed;
  }
  else {
    // Trigger a full chunk lookup
    // TODO: have a search tree here with logarithmic instead of linear complexity
    std::size_t chunk_ind = 0;
    for(chunk_ind = 0; chunk_ind < m_chunk_index.size(); chunk_ind++) {
      if(chunkContainsInds(m_chunk_index[chunk_ind], inds)) {
	m_chunk_last_accessed = chunk_ind;
	return chunk_ind;
      }
    }   
  }

  std::cout << "HHHHHHH" << std::endl;
  std::cout << "No chunk for index:" << std::endl;
  inds.print();
  std::cout << "HHHHHHH" << std::endl;
  
  throw ChunkNotFoundError();
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
      m_chunk_cache.insert({chunk_ind, dense_t::FromSparseFile(ifs)});
      // m_chunk_cache.insert({chunk_ind, dense_t::From(m_ser.template deserialize<sparse_t>(ifs))});
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
