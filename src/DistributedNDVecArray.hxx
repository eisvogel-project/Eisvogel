#include <uuid/uuid.h>
#include "DistributedNDVecArray.hh"

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
CacheEntry<ArrayT, T, dims, vec_dims>& CacheEntry<ArrayT, T, dims, vec_dims>::operator=(std::tuple<const metadata_t&, const chunk_t&, const status_t&> other) {

  chunk_meta = std::get<const metadata_t&>(other);
  chunk_data = std::get<const chunk_t&>(other);  
  op_to_perform = std::get<const status_t&>(other); 
  
  return *this;
}

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::ChunkCache(std::filesystem::path workdir, std::size_t cache_size, std::size_t init_cache_el_linear_size,
						  std::size_t initial_buffer_size) :
  ChunkCache(workdir, cache_size, Vector<std::size_t, dims>(init_cache_el_linear_size),
	     Vector<std::size_t, dims>(stor::INFTY), initial_buffer_size) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::ChunkCache(std::filesystem::path workdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
						  const Vector<std::size_t, dims>& streamer_chunk_size, std::size_t initial_buffer_size) :
  m_workdir(std::filesystem::absolute(workdir)), m_cache(cache_size, init_cache_el_shape), m_streamer(initial_buffer_size), m_streamer_chunk_size(streamer_chunk_size) {

  // create directory if it does not yet exist
  if(!std::filesystem::exists(m_workdir)) {
    std::filesystem::create_directory(m_workdir);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkCache<ArrayT, T, dims, vec_dims>::~ChunkCache() {
  FlushCache();
}

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
      [[unlikely]];
      
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
    if(std::get<CacheStatus::Append>(status).axis != axis) {
      
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
void ChunkCache<ArrayT, T, dims, vec_dims>::FlushCache() {

  // Evict and descope all cache elements
  for(id_t cur_index : m_cache.contained_elements()) {
    cache_entry_t& entry = m_cache.evict(cur_index);
    descope_cache_element(entry);
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

  // TODO: add check here whether ChunkType is specified or all_null
  
  // Directly deserialize into the cache element
  std::filesystem::path chunk_path = get_abs_path(chunk_meta.filepath);
  assert(std::filesystem::exists(chunk_path));
  
  std::fstream ifs;
  ifs.open(chunk_path, std::ios::in | std::ios::binary);
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

  std::filesystem::path chunk_path = get_abs_path(cache_entry.chunk_meta.filepath);
  assert(std::filesystem::exists(chunk_path));
  
  std::fstream iofs;
  iofs.open(chunk_path, std::ios::in | std::ios::out | std::ios::binary);
  iofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);

  // TODO: add check here whether ChunkType is specified or all_null
  
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

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::descope_cache_element(cache_entry_t& cache_entry) {

  std::filesystem::path chunk_path = get_abs_path(cache_entry.chunk_meta.filepath);
  
  // Handle any outstanding operations before this cache element goes out of scope and may be overwritten
  status_t& status = cache_entry.op_to_perform;
  if(std::holds_alternative<CacheStatus::Serialize>(status)) {
    
    // Serialize this chunk
    std::fstream ofs;
    ofs.open(chunk_path, std::ios::out | std::ios::binary);
    ofs.seekp(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    m_streamer.serialize(ofs, cache_entry.chunk_data, m_streamer_chunk_size, stor::StreamerMode::null_suppressed);
    ofs.close();    
  }
  else if(std::holds_alternative<CacheStatus::Append>(status)) {

    assert(std::filesystem::exists(chunk_path));
    
    // Need to append this to the on-disk chunk
    std::fstream iofs;
    iofs.open(chunk_path, std::ios::in | std::ios::out | std::ios::binary);
    iofs.seekg(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    iofs.seekp(cache_entry.chunk_meta.start_pos, std::ios_base::beg);
    m_streamer.append_slice(iofs, cache_entry.chunk_data, std::get<CacheStatus::Append>(status).axis,
			    stor::StreamerMode::null_suppressed);
    iofs.close();
  }
  else if(std::holds_alternative<CacheStatus::Nothing>(status)) {
    // Nothing needs to be done here
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
std::filesystem::path ChunkCache<ArrayT, T, dims, vec_dims>::get_abs_path(const std::filesystem::path& chunk_path) {
  return m_workdir / chunk_path;
}

// -------------

template <std::size_t dims>
ChunkMetadata<dims>::ChunkMetadata() :
  chunk_type(ChunkType::specified), filepath(""), chunk_id(0), start_pos(0), start_ind(0), end_ind(0), shape(0) { }

template <std::size_t dims>
ChunkMetadata<dims>::ChunkMetadata(const ChunkType& chunk_type, const std::filesystem::path& filepath, const id_t& chunk_id,
				   const Vector<std::size_t, dims> start_ind, const Vector<std::size_t, dims> shape) :
  chunk_type(chunk_type), filepath(filepath), chunk_id(chunk_id), start_ind(start_ind), shape(shape),
  end_ind(start_ind + shape) { }

template <std::size_t dims>
template <std::size_t axis>
void ChunkMetadata<dims>::GrowChunk(std::size_t shape_growth) {
  end_ind[axis] += shape_growth;
  shape[axis] += shape_growth;
}

template <std::size_t dims>
std::ostream& operator<<(std::ostream& stream, const ChunkMetadata<dims>& meta) {

  // TODO: make the output format nicer
  std::cout << "---------------------------------------\n";
  std::cout << " path:\t\t"      << meta.filepath  << "\n";
  std::cout << " chunk_id:\t\t"  << meta.chunk_id  << "\n";
  std::cout << " start_ind:\t\t" << meta.start_ind << "\n";
  std::cout << " end_ind:\t\t"   << meta.end_ind   << "\n";
  std::cout << " shape:\t\t"     << meta.shape     << "\n";
  std::cout << "---------------------------------------\n";
  
  return stream;
}

namespace stor{

  template <std::size_t dims>
  struct Traits<ChunkMetadata<dims>> {
    using type = ChunkMetadata<dims>;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<std::size_t>::serialize(stream, static_cast<std::size_t>(val.chunk_type));
      Traits<std::string>::serialize(stream, val.filepath.string());
      Traits<id_t>::serialize(stream, val.chunk_id);
      Traits<Vector<std::size_t, dims>>::serialize(stream, val.start_ind);
      Traits<Vector<std::size_t, dims>>::serialize(stream, val.shape);
    }

    static type deserialize(std::iostream& stream) {
      ChunkType chunk_type{Traits<std::size_t>::deserialize(stream)};
      std::filesystem::path filepath(Traits<std::string>::deserialize(stream));
      id_t chunk_id = Traits<id_t>::deserialize(stream);
      Vector<std::size_t, dims> start_ind = Traits<Vector<std::size_t, dims>>::deserialize(stream);
      Vector<std::size_t, dims> shape = Traits<Vector<std::size_t, dims>>::deserialize(stream);
      return ChunkMetadata<dims>(chunk_type, filepath, chunk_id, start_ind, shape);
    }
  };
}

// -------------

template <std::size_t dims>
ChunkIndex<dims>::ChunkIndex(std::filesystem::path index_path) : m_index_path(index_path), m_next_chunk_id(0), m_last_accessed_ind(0), m_shape(0) {

  // Already have an index file on disk, load it
  if(std::filesystem::exists(m_index_path)) {
    load_and_rebuild_index();
  }
}

template <std::size_t dims>
ChunkIndex<dims>::~ChunkIndex() {
  FlushIndex();
}

template <std::size_t dims>
ChunkIndex<dims>::metadata_t& ChunkIndex<dims>::RegisterChunk(const Vector<std::size_t, dims>& start_ind,
							      const Vector<std::size_t, dims>& shape) {

  uuid_t uuid_binary;
  uuid_generate_random(uuid_binary);
  char uuid_string[36];
  uuid_unparse(uuid_binary, uuid_string);
  std::filesystem::path filename(std::string(uuid_string) + ".bin");

  // build metadata object for this chunk
  id_t chunk_id = get_next_chunk_id();
  metadata_t chunk_meta(ChunkType::specified, filename, chunk_id, start_ind, shape);
  
  m_chunk_list.push_back(chunk_meta);
  return m_chunk_list.back();
}

template <std::size_t dims>
ChunkIndex<dims>::metadata_t& ChunkIndex<dims>::GetChunk(const Vector<std::size_t, dims>& ind) {

  assert(m_chunk_list.size() > 0);
    
  // First check if we're still in the most-recently read chunk  
  metadata_t& last_accessed_chunk = m_chunk_list[m_last_accessed_ind];
  if(is_in_chunk(last_accessed_chunk, ind)) {
    [[likely]];
    return last_accessed_chunk;
  }
  
  // If not, trigger full chunk lookup
  // TODO: replace this with lookup in the R-tree, which will be much more efficient
  for(std::size_t chunk_ind = 0; chunk_ind < m_chunk_list.size(); chunk_ind++) {
    metadata_t& cur_chunk = m_chunk_list[chunk_ind];
    if(is_in_chunk(cur_chunk, ind)) {
      m_last_accessed_ind = chunk_ind;
      return cur_chunk;
    }
  }

  throw std::logic_error("Error: no chunk found!");
}

template <std::size_t dims>
void ChunkIndex<dims>::FlushIndex() {

  // TODO: requires modifications after switch to R-tree
  std::fstream ofs;
  ofs.open(m_index_path, std::ios::out | std::ios::binary);
  stor::Traits<std::vector<metadata_t>>::serialize(ofs, m_chunk_list);
  ofs.close();
}

template <std::size_t dims>
ChunkIndex<dims>::shape_t ChunkIndex<dims>::GetShape() {

  // TODO: if this becomes a bottleneck, cache the shape and only invalidate it if an additional chunk gets added
  calculate_shape();
  return m_shape;
}

template <std::size_t dims>
void ChunkIndex<dims>::load_and_rebuild_index() {

  // TODO: requires modifications after switch to R-tree  
  m_chunk_list.clear();
  std::fstream ifs;
  ifs.open(m_index_path, std::ios::in | std::ios::binary);
  m_chunk_list = stor::Traits<std::vector<metadata_t>>::deserialize(ifs);
  ifs.close();
}

template <std::size_t dims>
bool ChunkIndex<dims>::is_in_chunk(const metadata_t& chunk, const Vector<std::size_t, dims>& ind) {
  return is_in_region(chunk.start_ind, chunk.shape, ind);
}

template <std::size_t dims>
bool ChunkIndex<dims>::is_in_region(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& shape,
				    const Vector<std::size_t, dims>& ind) {
  
  for(std::size_t i = 0; i < dims; i++) {
    
    // Efficient out-of-range comparison with a single branch
    if((ind[i] - start_ind[i]) >= shape[i]) {
      return false;
    }
  }  
  return true;
}

template <std::size_t dims>
id_t ChunkIndex<dims>::get_next_chunk_id() {
  return m_next_chunk_id++;
}

template <std::size_t dims>
void ChunkIndex<dims>::calculate_shape() {

  if(m_chunk_list.empty()) {
    return;
  }

  shape_t global_start_ind = get_start_ind();
  shape_t global_end_ind = get_end_ind();
  
  // check if the total inferred shape is consistent with the total number of elements contained in all chunks:
  // if so, then all chunks taken together define a contiguous region
  auto number_elements = [](const Vector<std::size_t, dims>& shape) -> std::size_t {
    return std::accumulate(shape.cbegin(), shape.cend(), 1, std::multiplies<std::size_t>());
  };

  std::size_t elements_from_shape = number_elements(global_end_ind - global_start_ind);
  
  std::size_t elements_in_chunks = 0;
  for(const metadata_t& cur_chunk : m_chunk_list) {
    elements_in_chunks += number_elements(cur_chunk.end_ind - cur_chunk.start_ind);
  }

  // have a contiguous region that is worth assigning a certain shape
  if(elements_in_chunks == elements_from_shape) {
    m_shape = global_end_ind - global_start_ind;
  }  
}

template <std::size_t dims>
Vector<std::size_t, dims> ChunkIndex<dims>::get_start_ind() {
  
  Vector<std::size_t, dims> start_ind(std::numeric_limits<std::size_t>::max());

  // TODO: this will be faster once we have the R-tree
  for(const metadata_t& cur_chunk : m_chunk_list) {
    for(std::size_t i = 0; i < dims; i++) {
      start_ind[i] = std::min(start_ind[i], cur_chunk.start_ind[i]);
    }
  }

  return start_ind;
}

template <std::size_t dims>
Vector<std::size_t, dims> ChunkIndex<dims>::get_end_ind() {
  
  Vector<std::size_t, dims> end_ind(std::numeric_limits<std::size_t>::min());

  // TODO: this will be faster once we have the R-tree
  for(const metadata_t& cur_chunk : m_chunk_list) {
    for(std::size_t i = 0; i < dims; i++) {
      end_ind[i] = std::max(end_ind[i], cur_chunk.end_ind[i]);
    }
  }
  
  return end_ind;
}

// -------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkLibrary<ArrayT, T, dims, vec_dims>::ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size,
						      std::size_t init_cache_el_linear_size) :
  ChunkLibrary(libdir, cache_size, Vector<std::size_t, dims>(init_cache_el_linear_size),
	       Vector<std::size_t, dims>(stor::INFTY)) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkLibrary<ArrayT, T, dims, vec_dims>::ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size,
						      const chunk_shape_t& init_cache_el_shape,
						      const Vector<std::size_t, dims>& streamer_chunk_size) :
  m_libdir(libdir), m_index_path(libdir / "index.bin"), m_index(m_index_path),
  m_cache(libdir, cache_size, init_cache_el_shape, streamer_chunk_size) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::RegisterChunk(const ind_t& start_ind, const chunk_t& chunk) {

  // Insert new metadata entry into chunk index
  metadata_t& meta = m_index.RegisterChunk(start_ind, chunk.GetShape());

  // Insert new chunk into the cache
  m_cache.RegisterNewChunk(meta, chunk);  
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::AppendSlice(const ind_t& start_ind, const chunk_t& slice) {
  
  // This is the index of an element in the (existing) chunk the slice should be appended to
  ind_t ind_existing_chunk = start_ind;
  ind_existing_chunk[axis] -= 1;
  
  // Get the metadata describing that chunk
  metadata_t& meta = m_index.GetChunk(ind_existing_chunk);
  
  // Perform the appending operation
  m_cache.template AppendSlice<axis>(meta, slice);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkLibrary<ArrayT, T, dims, vec_dims>::view_t ChunkLibrary<ArrayT, T, dims, vec_dims>::operator[](const ind_t& ind) {

  // Find chunk that holds the element with the required index
  metadata_t& meta = m_index.GetChunk(ind);
  
  // Retrieve the chunk
  const chunk_t& chunk = m_cache.RetrieveChunk(meta);

  // Fetch the element from the chunk
  return chunk[ind - meta.start_ind];
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::CalculateShape() {
  m_index.CalculateShape();
}

// ------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
DistributedNDVecArray<ArrayT, T, dims, vec_dims>::DistributedNDVecArray(std::filesystem::path workdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
									const Vector<std::size_t, dims>& streamer_chunk_size) :
  m_library(workdir, cache_size, init_cache_el_shape, streamer_chunk_size) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
DistributedNDVecArray<ArrayT, T, dims, vec_dims>::DistributedNDVecArray(std::filesystem::path workdir, std::size_t cache_size, std::size_t init_cache_el_linear_size) :
  m_library(workdir, cache_size, init_cache_el_linear_size) { }

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RegisterChunk(const ind_t& start_ind, const chunk_t& chunk) {
  m_library.RegisterChunk(start_ind, chunk);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
DistributedNDVecArray<ArrayT, T, dims, vec_dims>::view_t DistributedNDVecArray<ArrayT, T, dims, vec_dims>::operator[](const ind_t& ind) {
  return m_library[ind];
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::AppendSlice(const ind_t& start_ind, const chunk_t& slice) {
  m_library.template AppendSlice<axis>(start_ind, slice);
}
