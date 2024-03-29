#include "DistributedNDVecArray.hh"

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
CacheEntry<ArrayT, T, dims, vec_dims>& CacheEntry<ArrayT, T, dims, vec_dims>::operator=(std::tuple<metadata_t&, chunk_t&, CacheStatus&> other) {

  chunk_meta = std::get<metadata_t&>(other);
  chunk_data = std::get<chunk_t&>(other);  
  op_to_perform = std::get<CacheStatus&>(other); 
  
  return *this;
}

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::ChunkCache(std::size_t cache_size, const chunk_shape_t& init_cache_el_shape) :
  m_cache(cache_size, init_cache_el_shape, T()), m_streamer() { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::RegisterNewChunk(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data) {

  // Register a new chunk in the cache
  insert_into_cache(chunk_meta, chunk_data, CacheStatus::serialize);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
const ChunkCache<ArrayT, T, dims, vec_dims>::chunk_t& ChunkCache<ArrayT, T, dims, vec_dims>::RetrieveChunk(const chunk_meta_t& chunk_meta) {

  id_t index = chunk_meta.chunk_id;
  if(m_cache.contains(index)) {
    [[likely]]; // The point about caching is that it only makes sense when the access pattern is such that the same element is repeatedly needed many times

    // Requested chunk is contained in the cache, check if it needs to be brought up-to-date before returning a reference to it
    cache_entry_t& cached_chunk = m_cache.get(index);
    if((cached_chunk.op_to_perform != CacheStatus::nothing) &&
       (cached_chunk.op_to_perform != CacheStatus::serialize)) {
      
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
    insert_into_cache(chunk_meta, slice, CacheStatus::append);
  }
  else {

    // Cache already contains the chunk, how to proceed depends on its status
    cache_entry_t& cached_chunk = m_cache.get(index);

    switch(cached_chunk.op_to_perform) {
      
    case CacheStatus::append:
      // -> if the cache entry also has `append` as status, update the metadata and perform the concatenation in the cache, keep `append` as status
      break;

    case CacheStatus::serialize:
      // -> if the cache entry has `serialize` as status, update the metadata and perform the concatenation in the cache, keep `serialize` as status
      break;

    case CacheStatus::nothing:
      // -> if the cache entry has `nothing` as status, `free` it (so that it does not trigger any further lookups), remove its entry in the cache slot mapping
      //    and proceed as in 2)
      break;
    }
  } 
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::cache_entry_t& ChunkCache<ArrayT, T, dims, vec_dims>::deserialize_into_cache(const chunk_meta_t& chunk_meta) {

  if(!m_cache.has_free_slot()) {
    // Cache is full, evict oldest element and handle any outstanding operations
    cache_entry_t& oldest_entry = m_cache.evict_oldest_from_full_cache();
    oldest_entry.descope_entry();
  }
    
  // Now have free slot in the cache; fill it
  id_t index = chunk_meta.chunk_id;
  cache_entry_t& insert_location = m_cache.insert_ref_no_overwrite(index);

  insert_location.chunk_meta = chunk_meta;
  insert_location.op_to_perform = CacheStatus::nothing;  // this chunk is freshly read into the cache, nothing left to be done when it goes out of scope

  // Directly deserialize into the cache element
  assert(std::filesystem::exists(chunk_meta.filepath));
  
  std::fstream ifs;
  ifs.open(chunk_meta.filepath, std::ios::in | std::ios::binary);
  m_streamer.deserialize(ifs, insert_location.chunk_data);
  ifs.close();

  return insert_location;
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::insert_into_cache(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data, const CacheStatus& stat) {

  if(!m_cache.has_free_slot()) {
    // Cache is full, evict oldest element and handle any outstanding operations
    cache_entry_t& oldest_entry = m_cache.evict_oldest_from_full_cache();
    oldest_entry.descope_entry();
  }

  // Now have free slot in the cache, insert new element
  id_t index = chunk_meta.chunk_id;
  m_cache.insert_no_overwrite(index, std::forward_as_tuple(chunk_meta, chunk_data, stat));
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::sync_cache_element_for_read(cache_entry_t& cache_entry) {

  if(cache_entry.op_to_perform == CacheStatus::append) {
    // 3) `append`: perform appending to disk using the streamer, then reread into this slot (using the reference obtained at step 0)
  }
}

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void CacheEntry<ArrayT, T, dims, vec_dims>::descope_entry() {

  // 1) it becomes the oldest element in the cache and goes out of scope: determines what needs to be done to bring the 
  //    cache in sync with the on-disk representation of this chunk
  //    -> this is cleanly handled in the `descoper`
  //    *) nothing: can simply let it go out of scope,
  //    *) serialize: create the file listed in the metadata and serialize the cache entry into it, together with the metadata
  //    *) append: the file listed in the metadata already exists, append the chunk data and update the file metadata
  // 2) move status to `nothing` so that nothing happens to it when it is selected for another descope
  
  std::cout << "descoping cache entry" << std::endl;  
}

// -------------
