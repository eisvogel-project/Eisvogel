#include "DistributedNDVecArray.hh"

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
CacheEntry<ArrayT, T, dims, vec_dims>& CacheEntry<ArrayT, T, dims, vec_dims>::operator=(std::tuple<metadata_t&, chunk_t&, status_t&> other) {

  chunk_meta = std::get<metadata_t&>(other);
  chunk_data = std::get<chunk_t&>(other);  
  op_to_perform = std::get<status_t&>(other); 
  
  return *this;
}

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::ChunkCache(std::size_t cache_size, const chunk_shape_t& init_cache_el_shape, const Vector<std::size_t, dims>& streamer_chunk_size,
						  std::size_t initial_buffer_size) :
  m_cache(cache_size, init_cache_el_shape, T()), m_streamer(initial_buffer_size), m_streamer_chunk_size(streamer_chunk_size) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::RegisterNewChunk(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data) {

  // Register a new chunk in the cache
  insert_into_cache(chunk_meta, chunk_data, CacheStatus::Serialize());
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
const ChunkCache<ArrayT, T, dims, vec_dims>::chunk_t& ChunkCache<ArrayT, T, dims, vec_dims>::RetrieveChunk(const chunk_meta_t& chunk_meta) {

  id_t index = chunk_meta.chunk_id;
  if(m_cache.contains(index)) {
    [[likely]]; // The point about caching is that it only makes sense when the access pattern is such that the same element is repeatedly needed many times

    // Requested chunk is contained in the cache, check if it needs to be brought up-to-date before returning a reference to it
    cache_entry_t& cached_chunk = m_cache.get(index);
    if((!std::holds_alternative<CacheStatus::Nothing>(cached_chunk.op_to_perform)) &&
       (!std::holds_alternative<CacheStatus::Serialize>(cached_chunk.op_to_perform))) {
      
      // The element contained in the cache is not up-to-date; need to synchronize first
      sync_cache_element_for_read(cached_chunk);
    }
    
    // Requested element is contained in the cache and is now complete, directly return reference to its data
    return cached_chunk.chunk_data;
  }

  // Requested element is not yet contained in the cache, need to deserialize it first
  return deserialize_into_cache(chunk_meta).chunk_data;
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void ChunkCache<ArrayT, T, dims, vec_dims>::AppendSlice(chunk_meta_t& chunk_meta, const chunk_t& slice) {

  // Check if the concatenation along `axis` is compatible with the shape of the existing chunk and the slice
  if(!ArrayT<T, dims, vec_dims>::template ShapeAllowsConcatenation<axis>(chunk_meta.shape, slice.GetShape())) {
    throw std::runtime_error("Error: dimensions not compatible for concatenation!");
  }

  // Update the metadata
  std::size_t shape_growth = slice.GetShape()[axis];
  chunk_meta.template GrowChunk<axis>(shape_growth);
  
  id_t index = chunk_meta.chunk_id;
  if(!m_cache.contains(index)) {
    
    // Cache does not have the chunk to which the new slice should be appended; simply insert the slice as a new element into the cache using the `append` status
    // so that it will be appended to disk whenever it goes out of scope
    insert_into_cache(chunk_meta, slice, CacheStatus::Append(axis));
    return;
  }

  // Cache already contains the chunk, how to proceed depends on its status
  cache_entry_t& cached_chunk = m_cache.get(index);
  status_t& status = cached_chunk.op_to_perform;
  if(std::holds_alternative<CacheStatus::Nothing>(status)) {
    
    // The chunk is in the cache, but is marked as being in sync with the information on disk.
    // Remove it from the cache ...
    m_cache.evict(index);
    
    // ... and insert it with `append` status, as done above
    insert_into_cache(chunk_meta, slice, CacheStatus::Append(axis));
  }
  else if(std::holds_alternative<CacheStatus::Serialize>(status)) {
    
    // The cache already contains a chunk that is to be serialized from scratch; just concatenate in memory and write to disk whenever
    // this chunk is evicted from the cache
    cached_chunk.chunk_meta = chunk_meta;
    cached_chunk.chunk_data.template Append<axis>(slice);      
  }
  else if(std::holds_alternative<CacheStatus::Append>(status)) {
    
    // The cache already contains a chunk scheduled for concatenation with the on-disk chunk
    if(status.axis != axis) {
      
      // Concatenation axis changed, need to synchronize
      sync_cache_element_for_read(cached_chunk);
      cached_chunk.op_to_perform = CacheStatus::Append(axis); // record the new concatenation axis
    }

    // Now have the fully up-to date chunk in the cache, can append
    cached_chunk.chunk_meta = chunk_meta;
    cached_chunk.chunk_data.template Append<axis>(slice);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::cache_entry_t& ChunkCache<ArrayT, T, dims, vec_dims>::deserialize_into_cache(const chunk_meta_t& chunk_meta) {

  if(!m_cache.has_free_slot()) {
    // Cache is full, evict oldest element and handle any outstanding operations
    cache_entry_t& oldest_entry = m_cache.evict_oldest_from_full_cache();
    descope_cache_element(oldest_entry);
  }
    
  // Now have free slot in the cache; fill it
  id_t index = chunk_meta.chunk_id;
  cache_entry_t& insert_location = m_cache.insert_ref_no_overwrite(index);

  insert_location.chunk_meta = chunk_meta;
  insert_location.op_to_perform = CacheStatus::Nothing();  // this chunk is freshly read into the cache, nothing left to be done when it goes out of scope

  // Directly deserialize into the cache element
  assert(std::filesystem::exists(chunk_meta.filepath));
  
  std::fstream ifs;
  ifs.open(chunk_meta.filepath, std::ios::in | std::ios::binary);
  ifs.seekg(chunk_meta.start_pos, std::ios_base::beg);
  m_streamer.deserialize(ifs, insert_location.chunk_data);
  ifs.close();

  return insert_location;
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::insert_into_cache(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data, const status_t& stat) {

  if(!m_cache.has_free_slot()) {
    // Cache is full, evict oldest element and handle any outstanding operations
    cache_entry_t& oldest_entry = m_cache.evict_oldest_from_full_cache();
    descope_cache_element(oldest_entry);
  }

  // Now have free slot in the cache, insert new element
  id_t index = chunk_meta.chunk_id;
  m_cache.insert_no_overwrite(index, std::forward_as_tuple(chunk_meta, chunk_data, stat));
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::sync_cache_element_for_read(cache_entry_t& cache_entry) {

  std::filesystem::path chunk_path = cache_entry.chunk_meta.filepath;
  assert(std::filesystem::exists(chunk_path));
  
  std::fstream iofs;
  iofs.open(chunk_path, std::ios::in | std::ios::out | std::ios::binary);
  iofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);

  status_t& status = cache_entry.op_to_perform;
  if(std::holds_alternative<CacheStatus::Append>(status)) {
    
    // Append the cached chunk to the on-disk portion ...
    m_streamer.append_slice(iofs, cache_entry.chunk_data, std::get<CacheStatus::Append>(status).axis,
			    stor::StreamerMode::null_suppressed);
    
    // ... then rewind and deserialize the full chunk into the same cache location
    iofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    m_streamer.deserialize(iofs, cache_entry.chunk_data);
  }

  iofs.close();
}

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::descope_cache_element(cache_entry_t& cache_entry) {

  std::filesystem::path chunk_path = cache_entry.chunk_meta.filepath;
  assert(std::filesystem::exists(chunk_path));
  
  // Handle any outstanding operations before this cache element goes out of scope and may be overwritten
  status_t& status = cache_entry.op_to_perform;
  if(std::holds_alternative<CacheStatus::Serialize>(status)) {
    
    // Serialize this chunk
    std::fstream ofs;
    ofs.open(chunk_path, std::ios::out | std::ios::binary);
    ofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    m_streamer.serialize(ofs, cache_entry.chunk_data, m_streamer_chunk_size, stor::StreamerMode::null_suppressed);
    ofs.close();    
  }
  else if(std::holds_alternative<CacheStatus::Append>(status)) {
    
    // Need to append this to the on-disk chunk
    std::fstream iofs;
    iofs.open(chunk_path, std::ios::in | std::ios::out | std::ios::binary);
    iofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    m_streamer.append_slice(iofs, cache_entry.chunk_data, std::get<CacheStatus::Append>(status).axis,
			    stor::StreamerMode::null_suppressed);
    iofs.close();
  }
  else if(std::holds_alternative<CacheStatus::Nothing>(status)) {
    // Nothing needs to be done here
  }
}

// -------------
