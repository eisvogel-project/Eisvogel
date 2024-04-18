#include <uuid/uuid.h>
#include <algorithm>
#include <utility>
#include "DistributedNDVecArray.hh"
#include "Eisvogel/IteratorUtils.hh"

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
  m_workdir(std::filesystem::absolute(workdir)), m_streamer_chunk_size(streamer_chunk_size), m_cache(cache_size, init_cache_el_shape),
  m_streamer(initial_buffer_size) { }

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
void ChunkCache<ArrayT, T, dims, vec_dims>::RemoveChunk(const chunk_meta_t& chunk_meta) {

  // Remove the element from the cache whatever its current status is
  id_t index = chunk_meta.chunk_id;
  if(m_cache.contains(index)) {
    m_cache.evict(index);
  }

  // Also remove the chunk from disk
  if(std::filesystem::exists(chunk_meta.filepath)) {
    std::filesystem::remove(chunk_meta.filepath);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::ReplaceChunk(const chunk_meta_t& chunk_meta, const chunk_meta_t& new_chunk_meta, const chunk_t& new_chunk_data) {

  // replacing a chunk means that the chunk index does not change, but only other metadata or chunk data does
  assert(chunk_meta.chunk_id == new_chunk_meta.chunk_id);
  
  RemoveChunk(chunk_meta);
  RegisterNewChunk(new_chunk_meta, new_chunk_data);
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
      sync_cache_element_with_disk(cached_chunk);
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
      sync_cache_element_with_disk(cached_chunk);
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
void ChunkCache<ArrayT, T, dims, vec_dims>::MoveCache(std::filesystem::path new_workdir) {

  // Shut down the cache and make sure everything is safely on disk
  FlushCache();

  // Prepare the new working directory if it does not exist already
  if(!std::filesystem::exists(new_workdir)) {
    std::filesystem::create_directory(new_workdir);
  }
  
  // Move all the files to the new location
  for(auto const& cur_entry : std::filesystem::directory_iterator(m_workdir)) {
    if(cur_entry.path().extension() != m_suffix) {
      continue;
    }
    std::filesystem::copy(cur_entry, new_workdir);
    std::filesystem::remove(cur_entry);
  }

  // Update the working directory
  m_workdir = new_workdir;
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::ClearCache() {

  // Shut down the cache
  FlushCache();

  // Delete all files
  for(auto const& cur_entry : std::filesystem::directory_iterator(m_workdir)) {
    if(cur_entry.path().extension() != m_suffix) {
      continue;
    }
    std::filesystem::remove(cur_entry);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::ImportCache(std::filesystem::path workdir) {

  // Simply copy all the chunk files from the other work directory
  for(auto const& cur_entry : std::filesystem::directory_iterator(workdir)) {
    if(cur_entry.path().extension() != m_suffix) {
      continue;
    }
    std::filesystem::copy(cur_entry, m_workdir);
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

  // cache must not already contain this element
  id_t index = chunk_meta.chunk_id;
  assert(!m_cache.contains(index));
  
  if(!m_cache.has_free_slot()) {
    // Cache is full, evict oldest element and handle any outstanding operations
    cache_entry_t& oldest_entry = m_cache.evict_oldest_from_full_cache();
    descope_cache_element(oldest_entry);
  }

  // Now have free slot in the cache, insert new element
  m_cache.insert_no_overwrite(index, std::forward_as_tuple(chunk_meta, chunk_data, stat));
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkCache<ArrayT, T, dims, vec_dims>::sync_cache_element_with_disk(cache_entry_t& cache_entry) {

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

  cache_entry.op_to_perform = CacheStatus::Nothing();
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
  std::filesystem::path full_chunk_path = chunk_path;
  full_chunk_path += m_suffix;
  return m_workdir / full_chunk_path;
}

// -------------

template <std::size_t dims>
ChunkMetadata<dims>::ChunkMetadata() :
  chunk_type(ChunkType::specified), filepath(""), chunk_id(0), start_pos(0), shape(0), start_ind(0), end_ind(0), overlap(0), loc_ind_offset(0) { }

template <std::size_t dims>
ChunkMetadata<dims>::ChunkMetadata(const ChunkType& chunk_type, const std::filesystem::path& filepath, const id_t& chunk_id,
				   const Vector<std::size_t, dims> start_ind, const Vector<std::size_t, dims> shape,
				   std::size_t overlap) :
  chunk_type(chunk_type), filepath(filepath), chunk_id(chunk_id), start_pos(0), shape(shape),
  start_ind(start_ind), end_ind(start_ind + shape), overlap(overlap), loc_ind_offset(start_ind - overlap) {

  // Note: unsigned-integer underflow in `loc_ind_offset` can happen, but is no problem here:
  // when a chunk element is accessed, the index will be `global_index - loc_ind_offset`.
}

template <std::size_t dims>
template <std::size_t axis>
void ChunkMetadata<dims>::GrowChunk(std::size_t shape_growth) {
  end_ind[axis] += shape_growth;
  shape[axis] += shape_growth;
}

template <std::size_t dims>
template <std::size_t axis_1, std::size_t axis_2>
void ChunkMetadata<dims>::SwapAxes() {
  std::swap(start_ind[axis_1], start_ind[axis_2]);
  std::swap(end_ind[axis_1], end_ind[axis_2]);
  std::swap(shape[axis_1], shape[axis_2]);
  std::swap(loc_ind_offset[axis_1], loc_ind_offset[axis_2]);
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
  std::cout << " loc_ind_offset:\t\t"   << meta.loc_ind_offset   << "\n";
  std::cout << " overlap:\t\t" << meta.overlap << "\n";
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
      Traits<std::size_t>::serialize(stream, val.overlap);
    }

    static type deserialize(std::iostream& stream) {
      ChunkType chunk_type{Traits<std::size_t>::deserialize(stream)};
      std::filesystem::path filepath(Traits<std::string>::deserialize(stream));
      id_t chunk_id = Traits<id_t>::deserialize(stream);
      Vector<std::size_t, dims> start_ind = Traits<Vector<std::size_t, dims>>::deserialize(stream);
      Vector<std::size_t, dims> shape = Traits<Vector<std::size_t, dims>>::deserialize(stream);
      std::size_t overlap = Traits<std::size_t>::deserialize(stream);
      return ChunkMetadata<dims>(chunk_type, filepath, chunk_id, start_ind, shape, overlap);
    }
  };
}

// -------------

template <std::size_t dims>
ChunkIndex<dims>::ChunkIndex(std::filesystem::path index_path) :
  m_next_chunk_id(0), m_index_path(index_path), m_index_metadata_valid(false), m_shape(0), m_start_ind(0), m_end_ind(0), m_chunk_list(), m_last_accessed_ind(0) {

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
							      const Vector<std::size_t, dims>& shape, std::size_t overlap) {

  uuid_t uuid_binary;
  uuid_generate_random(uuid_binary);
  char uuid_string[36];
  uuid_unparse(uuid_binary, uuid_string);
  std::filesystem::path filename = std::string(uuid_string);

  // adding a new chunk invalidates the previous shape
  invalidate_cached_index_metadata();
  
  // build metadata object for this chunk
  id_t chunk_id = get_next_chunk_id();
  metadata_t chunk_meta(ChunkType::specified, filename, chunk_id, start_ind, shape, overlap);
  
  m_chunk_list.push_back(chunk_meta);
  return m_chunk_list.back();
}

template <std::size_t dims>
const ChunkIndex<dims>::metadata_t& ChunkIndex<dims>::GetChunk(const Vector<std::size_t, dims>& ind) const {
  return find_chunk_by_index(ind);
}

template <std::size_t dims>
ChunkIndex<dims>::metadata_t& ChunkIndex<dims>::GetChunk(const Vector<std::size_t, dims>& ind) {

  // since we are returning a non-const reference to the chunk metadata, assume that the caller indends to modify it and that we need to
  // recompute any cached quantities describing the full chunk index
  invalidate_cached_index_metadata(); 
  return find_chunk_by_index(ind);
}

template <std::size_t dims>
std::vector<std::reference_wrapper<const ChunkMetadata<dims>>> ChunkIndex<dims>::GetChunks(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind) {

  std::vector<std::reference_wrapper<const metadata_t>> chunks_found;
  
  // TODO: this will be significantly faster once we can do the lookup in the R-tree instead of through a linear search
  for(const metadata_t& cur_chunk : m_chunk_list) {
    if(chunk_overlaps_region(cur_chunk, start_ind, end_ind)) {
      chunks_found.push_back(cur_chunk);
    }
  }
  return chunks_found;
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
void ChunkIndex<dims>::MoveIndex(std::filesystem::path new_index_path) {

  assert(new_index_path != m_index_path);

  // remove the old index file on disk
  if(std::filesystem::exists(m_index_path)) {
    std::filesystem::remove(m_index_path);
  }
  
  // reset the path and flush the index to the new location
  m_index_path = new_index_path;  
  FlushIndex();
}

template <std::size_t dims>
void ChunkIndex<dims>::ClearIndex() {

  if(std::filesystem::exists(m_index_path)) {
    std::filesystem::remove(m_index_path);
  }

  // reset and clear everything
  invalidate_cached_index_metadata();
  m_last_accessed_ind = 0;
  m_shape = 0;
  m_chunk_list.clear();
}

template <std::size_t dims>
void ChunkIndex<dims>::ImportIndex(std::filesystem::path index_path) {

  invalidate_cached_index_metadata();

  // Load the new index entries
  std::vector<metadata_t> index_entries = load_index_entries(index_path);

  // Re-enumerate the newly imported entries
  for(metadata_t& cur_meta : index_entries) {
    cur_meta.chunk_id = get_next_chunk_id();    
  }
  
  m_chunk_list.insert(m_chunk_list.end(),
		      std::make_move_iterator(index_entries.begin()),
		      std::make_move_iterator(index_entries.end())
		      );

  // TODO: rebuild the R-tree based on the new list of chunks
}

template <std::size_t dims>
ChunkIndex<dims>::shape_t ChunkIndex<dims>::get_start_ind() {
  if(!m_index_metadata_valid) {
    calculate_and_cache_index_metadata();
  }  
  return m_start_ind;
}

template <std::size_t dims>
ChunkIndex<dims>::shape_t ChunkIndex<dims>::get_end_ind() {
  if(!m_index_metadata_valid) {
    calculate_and_cache_index_metadata();
  }  
  return m_end_ind;
}

template <std::size_t dims>
ChunkIndex<dims>::shape_t ChunkIndex<dims>::GetShape() {
  if(!m_index_metadata_valid) {
    calculate_and_cache_index_metadata();
  }  
  return m_shape;
}

template <std::size_t dims>
std::vector<ChunkMetadata<dims>> ChunkIndex<dims>::load_index_entries(std::filesystem::path index_path) {

  std::cout << "loading index entries from " << index_path << std::endl;
  
  std::fstream ifs;
  ifs.open(index_path, std::ios::in | std::ios::binary);
  std::vector<metadata_t> index_entries = stor::Traits<std::vector<metadata_t>>::deserialize(ifs);
  ifs.close();

  std::cout << "found " << index_entries.size() << " entries" << std::endl;
  
  return index_entries;
}

template <std::size_t dims>
void ChunkIndex<dims>::load_and_rebuild_index() {

  m_chunk_list.clear();
  m_chunk_list = load_index_entries(m_index_path);

  // TODO: rebuild the R-tree based on the new list of chunks
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
bool ChunkIndex<dims>::chunk_overlaps_region(const metadata_t& chunk, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind) {

  for(std::size_t i = 0; i < dims; i++) {

    // the `start_ind` of the specified region must lie "to the left of" the chunk endpoint in all directions
    if(start_ind[i] >= chunk.end_ind[i]) {
      return false;
    }

    // the `end_ind` of the specified region must lie "to the right of" the chunk endpoint in all directions
    if(end_ind[i] <= chunk.start_ind[i]) {
      return false;
    }
  }  
  return true;
}

template <std::size_t dims>
Vector<std::size_t, dims> ChunkIndex<dims>::get_overlap_start_ind(const metadata_t& chunk, const Vector<std::size_t, dims>& start_ind) {

  // find element-wise std::max between `start_ind` (belonging to the range) and `chunk.start_ind`
  Vector<std::size_t, dims> overlap_start_ind;
  std::transform(std::execution::unseq, start_ind.begin(), start_ind.end(), chunk.start_ind.begin(), overlap_start_ind.begin(),
		 [](auto a, auto b){return std::max(a, b);});
  
  return overlap_start_ind;
}

template <std::size_t dims>
Vector<std::size_t, dims> ChunkIndex<dims>::get_overlap_end_ind(const metadata_t& chunk, const Vector<std::size_t, dims>& end_ind) {
  
  // find element-wise std::min between `end_ind` (belonging to the range) and `chunk.end_ind`
  Vector<std::size_t, dims> overlap_end_ind;
  std::transform(std::execution::unseq, end_ind.begin(), end_ind.end(), chunk.end_ind.begin(), overlap_end_ind.begin(),
		 [](auto a, auto b){return std::min(a, b);});
  
  return overlap_end_ind;
}

template <std::size_t dims>
id_t ChunkIndex<dims>::get_next_chunk_id() {
  return m_next_chunk_id++;
}

template <std::size_t dims>
ChunkIndex<dims>::metadata_t& ChunkIndex<dims>::find_chunk_by_index(const Vector<std::size_t, dims>& ind) {

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
void ChunkIndex<dims>::invalidate_cached_index_metadata() {
  m_index_metadata_valid = false;
}

template <std::size_t dims>
void ChunkIndex<dims>::calculate_and_cache_index_metadata() {

  if(m_chunk_list.empty()) {
    m_start_ind = 0;
    m_end_ind = 0;
    m_shape = 0;
    m_index_metadata_valid = true;
    return;
  }

  // calculate start and end indices
  m_start_ind = calculate_start_ind();
  m_end_ind = calculate_end_ind();
  
  // check if the total inferred shape is consistent with the total number of elements contained in all chunks:
  // if so, then all chunks taken together define a contiguous region
  auto number_elements = [](const Vector<std::size_t, dims>& shape) -> std::size_t {
    return std::accumulate(shape.cbegin(), shape.cend(), 1, std::multiplies<std::size_t>());
  };

  std::size_t elements_from_shape = number_elements(m_end_ind - m_start_ind);
  
  std::size_t elements_in_chunks = 0;
  for(const metadata_t& cur_chunk : m_chunk_list) {
    elements_in_chunks += number_elements(cur_chunk.end_ind - cur_chunk.start_ind);
  }

  // have a contiguous region that is worth assigning a certain shape
  if(elements_in_chunks == elements_from_shape) {
    m_shape = m_end_ind - m_start_ind;
  }
  else {
    m_shape = 0;
  }

  m_index_metadata_valid = true;
}

template <std::size_t dims>
Vector<std::size_t, dims> ChunkIndex<dims>::calculate_start_ind() {
  
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
Vector<std::size_t, dims> ChunkIndex<dims>::calculate_end_ind() {
  
  Vector<std::size_t, dims> end_ind(std::numeric_limits<std::size_t>::min());

  // TODO: this will be faster once we have the R-tree
  for(const metadata_t& cur_chunk : m_chunk_list) {
    for(std::size_t i = 0; i < dims; i++) {
      end_ind[i] = std::max(end_ind[i], cur_chunk.end_ind[i]);
    }
  }
  
  return end_ind;
}

template <std::size_t dims>
auto ChunkIndex<dims>::begin() {
  // This returns a non-const iterator, make sure to invalidate any cached values
  invalidate_cached_index_metadata();
  return m_chunk_list.begin();
}

template <std::size_t dims>
auto ChunkIndex<dims>::end() {
  invalidate_cached_index_metadata();
  return m_chunk_list.end();
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
  m_libdir(libdir), m_index_path(libdir / m_index_filename), m_index(m_index_path),
  m_cache(libdir, cache_size, init_cache_el_shape, streamer_chunk_size) {

  // create directory if this does not exist already
  if(!std::filesystem::exists(libdir)) {
    std::filesystem::create_directory(libdir);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::RegisterChunk(const chunk_t& chunk, const ind_t& start_ind, const ind_t& end_ind, std::size_t overlap) {

  assert(chunk.GetShape() == end_ind - start_ind + 2 * overlap);
  
  // Insert new metadata entry into chunk index
  metadata_t& meta = m_index.RegisterChunk(start_ind, end_ind - start_ind, overlap);

  // Insert new chunk into the cache
  m_cache.RegisterNewChunk(meta, chunk);  
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::AppendSlice(const ind_t& start_ind, const chunk_t& slice) {

  std::cout << "BBBB appending, start_ind = " << start_ind << ", slice_shape = " << slice.GetShape() << std::endl;
  
  // This is the index of an element in the (existing) chunk the slice should be appended to
  ind_t ind_existing_chunk = start_ind;
  ind_existing_chunk[axis] -= 1;
  
  // Get the metadata describing that chunk
  metadata_t& meta = m_index.GetChunk(ind_existing_chunk);

  std::cout << "BBBB before append: \n" << meta << std::endl;
  
  // Perform the appending operation
  m_cache.template AppendSlice<axis>(meta, slice);

  std::cout << "BBBB after append: \n" << meta << std::endl;

  std::cout << "BBBB retrieved_shape = " << m_cache.RetrieveChunk(meta).GetShape() << ", meta.shape = " << meta.shape << std::endl;
  std::cout << "BBBB retrieved_shape = " << m_cache.RetrieveChunk(meta).GetShape() << ", meta.shape = " << meta.shape << std::endl;
  
  assert(meta.shape == m_cache.RetrieveChunk(meta).GetShape());
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
ChunkLibrary<ArrayT, T, dims, vec_dims>::view_t ChunkLibrary<ArrayT, T, dims, vec_dims>::operator[](const ind_t& ind) {

  // Find chunk that holds the element with the required index
  const metadata_t& meta = m_index.GetChunk(ind);
  
  // Retrieve the chunk
  const chunk_t& chunk = m_cache.RetrieveChunk(meta);

  // Convert to chunk-local index and fetch the element from the chunk
  return chunk[ind - meta.loc_ind_offset];
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::FillArray(chunk_t& array, const ind_t& input_start, const ind_t& input_end, const ind_t& output_start) {

  // Make sure the full range is available in the target array
  assert(array.has_index(output_start));
  assert(array.has_index(output_start + (input_end - input_start) - 1));
  
  // Go through all chunks, find overlaps with the targeted region, and fill these into the target array
  
  // Get chunks that (perhaps with only part of their elements) contribute to the targeted range
  std::vector<std::reference_wrapper<const metadata_t>> required_chunks = m_index.GetChunks(input_start, input_end);
  for(const metadata_t& chunk_meta : required_chunks) {

    const chunk_t& chunk = m_cache.RetrieveChunk(chunk_meta);

    // These are the indices where the current `chunk` overlaps with the specified range
    ind_t chunk_overlap_start_ind = ChunkIndex<dims>::get_overlap_start_ind(chunk_meta, input_start);
    ind_t chunk_overlap_end_ind = ChunkIndex<dims>::get_overlap_end_ind(chunk_meta, input_end);

    std::cout << chunk_meta << std::endl;
    std::cout << "chunk.shape = " << chunk.GetShape() << std::endl;
    std::cout << "chunk_overlap_start_ind = " << chunk_overlap_start_ind << std::endl;
    std::cout << "chunk_overlap_end_ind = " << chunk_overlap_end_ind << std::endl;
    std::cout << "input_start_chunk_local = " << chunk_overlap_start_ind - chunk_meta.loc_ind_offset << std::endl;
    
    // copy the overlapping range from the `chunk` into the destination `array`
    array.fill_from(chunk,
		    chunk_overlap_start_ind - chunk_meta.loc_ind_offset, // `input_start` in chunk-local indices
		    chunk_overlap_end_ind - chunk_meta.loc_ind_offset, // `input_end` in chunk-local indices
		    chunk_overlap_start_ind - input_start + output_start); // `output_start` in output-array-local indices
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis_1, std::size_t axis_2>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::SwapAxes() {

  auto swapper = [](metadata_t& new_chunk_meta, const chunk_t& chunk_data, chunk_t& new_chunk_data) {
    
    // swap axes in metadata entries
    new_chunk_meta.template SwapAxes<axis_1, axis_2>();
    
    // swap data in the chunk
    chunk_data.template SwapAxes<axis_1, axis_2>(new_chunk_data);
  };
  map_over_chunks(swapper);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::MoveLibrary(std::filesystem::path new_libdir) {

  // prepare the new directory
  if(!std::filesystem::exists(new_libdir)) {
    std::filesystem::create_directory(new_libdir);
  }
  
  // move everything to a new location
  m_cache.MoveCache(new_libdir);
  m_index.MoveIndex(new_libdir / m_index_filename);

  // delete the old directory
  std::filesystem::remove(m_libdir);

  m_libdir = new_libdir;
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::ClearLibrary() {

  // clear all contents, but do not delete the containing directory
  m_index.ClearIndex();
  m_cache.ClearCache();
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::DeleteLibrary() {

  // clear the library ...
  ClearLibrary();

  // ... and delete the directory
  std::filesystem::remove(m_libdir);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::FlushLibrary() {
  m_index.FlushIndex();
  m_cache.FlushCache();
}


template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::ImportLibrary(std::filesystem::path libdir) {

  // import and merge both the chunk index and the "cache"
  m_index.ImportIndex(libdir / m_index_filename);
  m_cache.ImportCache(libdir);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void ChunkLibrary<ArrayT, T, dims, vec_dims>::index_loop_over_elements(CallableT&& worker) {
  Vector<std::size_t, dims> start_ind(0);
  index_loop_over_elements(start_ind, m_index.GetShape(), worker);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void ChunkLibrary<ArrayT, T, dims, vec_dims>::index_loop_over_elements(const ind_t& start_ind, const ind_t& end_ind, CallableT&& worker) {

  // Get chunks that (perhaps with only part of their elements) contribute to the targeted range
  std::vector<std::reference_wrapper<const metadata_t>> required_chunks = m_index.GetChunks(start_ind, end_ind);
  for(const metadata_t& chunk_meta : required_chunks) {
    
    const chunk_t& chunk = m_cache.RetrieveChunk(chunk_meta);

    // These are the indices where the current `chunk` overlaps with the specified range
    ind_t chunk_overlap_start_ind = ChunkIndex<dims>::get_overlap_start_ind(chunk_meta, start_ind);
    ind_t chunk_overlap_end_ind = ChunkIndex<dims>::get_overlap_end_ind(chunk_meta, end_ind);
    
    auto chunk_worker = [&](const Vector<std::size_t, dims>& ind_within_chunk, const view_t& elem) {
      worker(chunk_meta.loc_ind_offset + ind_within_chunk, elem);
    };
    chunk.index_loop_over_elements(chunk_overlap_start_ind - chunk_meta.loc_ind_offset, chunk_overlap_end_ind - chunk_meta.loc_ind_offset, // chunk-local start and end indices
				   chunk_worker);
  }
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
void ChunkLibrary<ArrayT, T, dims, vec_dims>::map_over_chunks(CallableT&& worker) {

  chunk_t new_chunk_data(Vector<std::size_t, dims>(1));
  metadata_t new_chunk_meta;
    
  for(metadata_t& chunk_meta : m_index) {

    const chunk_t& chunk_data = m_cache.RetrieveChunk(chunk_meta);
    
    new_chunk_meta = chunk_meta; // initialize the new metadata with the current one

    // ask the worker to fill in the new chunk and modify the chunk metadata
    worker(new_chunk_meta, chunk_data, new_chunk_data);
    
    // use the new chunk to replace the old one in the cache
    m_cache.ReplaceChunk(chunk_meta, new_chunk_meta, new_chunk_data);

    // record the modification also in the chunk index
    chunk_meta = new_chunk_meta;
  }  
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
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RegisterChunk(const chunk_t& chunk, const ind_t& start_ind) {
  std::size_t overlap = 0;
  ind_t end_ind = start_ind + chunk.GetShape();
  m_library.RegisterChunk(chunk, start_ind, end_ind, overlap);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RegisterChunk(const chunk_t& chunk, const ind_t& start_ind, const ind_t& end_ind, std::size_t overlap) {
  m_library.RegisterChunk(chunk, start_ind, end_ind, overlap);
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

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
constexpr void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::index_loop_over_elements(CallableT&& worker) {
  m_library.index_loop_over_elements(worker);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::FillArray(chunk_t& array, const ind_t& start_ind, const ind_t& end_ind, const ind_t& output_start) {
  m_library.FillArray(array, start_ind, end_ind, output_start);
}


template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <std::size_t axis_1, std::size_t axis_2>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::SwapAxes() {
  m_library.template SwapAxes<axis_1, axis_2>();
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RebuildChunksPartial(const ind_t& start_ind, const ind_t& end_ind,
									    const ind_t& requested_chunk_shape, std::filesystem::path outdir) {

  // no overlap requested here
  auto error_on_evaluation = []() {
    throw std::logic_error("This should never be encountered!");
  };
  std::size_t overlap = 0;
  RebuildChunksPartial(start_ind, end_ind, requested_chunk_shape, outdir, overlap, error_on_evaluation);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class BoundaryCallableT>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RebuildChunksPartial(const ind_t& start_ind, const ind_t& end_ind,
									    const ind_t& requested_chunk_shape, std::filesystem::path outdir,
									    std::size_t overlap, BoundaryCallableT&& boundary_evaluator) {

  // prepare the output library containing the rebuilt chunks
  ind_t streamer_chunk_size(stor::INFTY);
  streamer_chunk_size[0] = 1;
  std::size_t cache_depth = 1; // no need for a large cache here, will just add one chunk at a time
  ChunkLibrary<ArrayT, T, dims, vec_dims> rebuilt_library(outdir, cache_depth, requested_chunk_shape, streamer_chunk_size);

  // buffer where the new chunks will be filled
  ArrayT<T, dims, vec_dims> chunk_buffer(requested_chunk_shape);

  // ------------------------------
  //  Note: std library containers like `std::array` strictly use unsigned types for indexing, and so do `NDVecArray` and `Vector`.
  //  Signed indices are useful for what follows, and so there will be a bit of back-and-forth type casting below, please bear with us.
  // ------------------------------

  ind_t global_start_ind = m_library.GetStartInd();
  ind_t global_end_ind = m_library.GetEndInd();
  
  auto rebuilder = [&](const ind_t& chunk_start_ind, const ind_t& chunk_end_ind) {

    // Compute the start and end indices of the chunks extended by the corresponding overlap ...
    // ... which means that these indices may be outside the bounds of the original array.    
    signed_ind_t extended_chunk_start_ind = chunk_start_ind.template as_type<int>() - overlap;
    signed_ind_t extended_chunk_end_ind = chunk_end_ind.template as_type<int>() + overlap;
    chunk_shape_t extended_chunk_shape = (chunk_end_ind - chunk_start_ind) + 2 * overlap;

    // Determine the index range that overlaps with the original array from which the contained values can be copied over
    ind_t fill_start_ind = VectorUtils::max(global_start_ind, extended_chunk_start_ind);
    ind_t fill_end_ind = VectorUtils::min(global_end_ind, extended_chunk_end_ind);

    // position in the output chunk where the copying should start
    ind_t output_start = (fill_start_ind.template as_type<int>() - extended_chunk_start_ind).template as_type<std::size_t>();
    
    // fill rebuilt chunk into local buffer ...    
    chunk_buffer.resize(extended_chunk_shape);    
    FillArray(chunk_buffer, fill_start_ind, fill_end_ind, output_start);
    
    // ... if introducing the overlap takes us beyond the global shape of the array, manually iterate over the boundary faces and ask the boundary evaluator
    // to fill in the field values there
    auto boundary_filler = [&](Vector<int, dims>& boundary_ind) {
      ind_t local_ind = (boundary_ind - extended_chunk_start_ind).template as_type<std::size_t>();
      boundary_evaluator(*this, std::forward<Vector<int, dims>>(boundary_ind), chunk_buffer[local_ind]);

      // std::cout << "set to ";
      // for(auto cur: chunk_buffer[local_ind]) {
      // 	std::cout << cur << " ";
      // }
      // std::cout << std::endl;
      
    };
    index_loop_over_penetrating_chunk_elements(global_start_ind, global_end_ind, extended_chunk_start_ind, extended_chunk_end_ind, boundary_filler);
    
    // ... and finally register the thus constructed chunk in the new library under the original start index, and pass on the information about the overlap
    rebuilt_library.RegisterChunk(chunk_buffer, chunk_start_ind, chunk_end_ind, overlap);
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, requested_chunk_shape, rebuilder);

  // Ensure that everything is written to disk before we finish here
  rebuilt_library.FlushLibrary();
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class CallableT>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::index_loop_over_penetrating_chunk_elements(const ind_t& start_ind, const ind_t& end_ind,
												  const signed_ind_t& chunk_start_ind, const signed_ind_t& chunk_end_ind,
												  CallableT&& worker) {
  // Determine the range where the chunk intersects the specified array range
  ind_t intersect_start_ind = VectorUtils::max(start_ind, chunk_start_ind);
  ind_t intersect_end_ind = VectorUtils::min(end_ind, chunk_end_ind);

  // Check penetrating slices along all directions
  for(std::size_t dim = 0; dim < dims; dim++) {

    // Chunk starts "to the left" of the array along this direction ...
    if(std::cmp_less(chunk_start_ind[dim], start_ind[dim])) {

      // ... need to iterate over the region that starts at the `chunk_start_ind` and extends up to the
      // `intersect_start_ind` along the direction of `dim`
      signed_ind_t boundary_slice_end_ind = chunk_end_ind;
      boundary_slice_end_ind[dim] = intersect_start_ind[dim];      
      IteratorUtils::index_loop_over_elements(chunk_start_ind, boundary_slice_end_ind, worker);
    }

    // Chunk ends "to the right" of the array along this direction ...
    if(std::cmp_greater(chunk_end_ind[dim], end_ind[dim])) {

      // ... need to iterate over the region that ends at the `chunk_end_ind` and extends backwards
      // along `dim` to where the intersection begins
      signed_ind_t boundary_slice_start_ind = chunk_start_ind;
      boundary_slice_start_ind[dim] = intersect_end_ind[dim];
      IteratorUtils::index_loop_over_elements(boundary_slice_start_ind, chunk_end_ind, worker);
    }
  }  
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::Import(std::filesystem::path dir) {
  m_library.ImportLibrary(dir);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::Move(std::filesystem::path dest) {
  m_library.MoveLibrary(dest);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
template <class BoundaryCallableT>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RebuildChunks(const ind_t& requested_chunk_shape, std::filesystem::path tmpdir, std::size_t overlap,
								     BoundaryCallableT&& boundary_evaluator) {
  ind_t global_start_ind(0);
  shape_t global_shape = m_library.GetShape();

  std::cout << "global_shape in rebuilding chunks = " << global_shape << std::endl;
  
  RebuildChunksPartial(global_start_ind, global_shape, requested_chunk_shape, tmpdir, overlap, boundary_evaluator);

  // Clear the contents of the original library ...
  m_library.ClearLibrary();

  // ... pull in the contents of the rebuilt one ...
  m_library.ImportLibrary(tmpdir);
  m_library.FlushLibrary();

  // ... and delete the temporary directory
  std::filesystem::remove_all(tmpdir);
}

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
void DistributedNDVecArray<ArrayT, T, dims, vec_dims>::RebuildChunks(const ind_t& requested_chunk_shape,
								     std::filesystem::path tmpdir) {
  // no overlap requested here
  auto error_on_evaluation = []() {
    throw std::logic_error("This should never be encountered!");
  };
  std::size_t overlap = 0;
  RebuildChunks(requested_chunk_shape, tmpdir, overlap, error_on_evaluation);
}
