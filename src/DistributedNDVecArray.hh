#pragma once

#include <filesystem>
#include <fstream>
#include <tuple>
#include <vector>
#include <variant>

#include "Cache.hh"
#include "Vector.hh"
#include "NDVecArrayStreamer.hh"

enum class ChunkType : std::size_t {
  specified = 0,
  all_null = 1
};

// forward declarations (needed for `operator<<` decleared below)
template <std::size_t dims>
struct ChunkMetadata;

template <std::size_t dims>
std::ostream& operator<<(std::ostream& stream, const ChunkMetadata<dims>& meta);

template <std::size_t dims>
struct ChunkMetadata {

  using id_t = std::size_t;

  // default constructor
  ChunkMetadata();
  ChunkMetadata(const ChunkType& chunk_type, const std::filesystem::path& filepath, const id_t& chunk_id,
		const Vector<std::size_t, dims> start_ind, const Vector<std::size_t, dims> shape);

  friend std::ostream& operator<< <dims> (std::ostream& stream, const ChunkMetadata& meta);
  
  // Accessor that grows the recorded shape of a chunk
  template <std::size_t axis>
  void GrowChunk(std::size_t shape_growth);
  
  // type of this chunk
  ChunkType chunk_type;
  
  // path to the file where this chunk lives
  std::filesystem::path filepath;

  // unique id: does not change when slices are appended
  id_t chunk_id;
  
  // offset within the file where this chunk lives
  std::streampos start_pos;

  // indices of the elements marking the start and end element in the chunk ...
  Vector<std::size_t, dims> start_ind;
  Vector<std::size_t, dims> end_ind;

  // ... as well as the (precalculated) shape
  Vector<std::size_t, dims> shape;
};

// -------------------------

// provides fast coordinate-indexed lookup of chunks
template <std::size_t dims>
class ChunkIndex {

public:

  using metadata_t = ChunkMetadata<dims>;
  using shape_t = Vector<std::size_t, dims>;
  
public:

  // build from stored index file
  ChunkIndex(std::filesystem::path index_path);
  ~ChunkIndex();

  metadata_t& RegisterChunk(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& shape);
  
  // get chunk that contains the element with index `ind`
  metadata_t& GetChunk(const Vector<std::size_t, dims>& ind);

  // get chunks that, taken together, cover the rectangular region between `start_ind` and `end_ind`
  std::vector<metadata_t> GetChunks(const Vector<std::size_t, dims>& start_ind,
				    const Vector<std::size_t, dims>& end_ind);  

  void FlushIndex();

  shape_t GetShape();

private:

  void calculate_shape();
  Vector<std::size_t, dims> get_start_ind();
  Vector<std::size_t, dims> get_end_ind();
  
  id_t get_next_chunk_id();

  void load_and_rebuild_index();
  
  bool is_in_chunk(const metadata_t& chunk, const Vector<std::size_t, dims>& ind);
  bool is_in_region(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& shape,
		    const Vector<std::size_t, dims>& ind);
  
private:

  id_t m_next_chunk_id;
  std::filesystem::path m_index_path;
  
  shape_t m_shape;
  
  // TODO: use R-tree to quickly find the chunks in this list
  std::vector<metadata_t> m_chunk_list;
  std::size_t m_last_accessed_ind;
};

// -------------------------

namespace CacheStatus {

  struct Nothing { };
  struct Serialize { };
  
  struct Append {
    Append(std::size_t axis) : axis(axis) { };
    std::size_t axis;
  };  
}

// The elements stored in the cache
template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
struct CacheEntry {  

  using metadata_t = ChunkMetadata<dims>;
  using chunk_t = ArrayT<T, dims, vec_dims>;
  using shape_t = typename chunk_t::shape_t;
  using status_t = std::variant<CacheStatus::Nothing, CacheStatus::Serialize, CacheStatus::Append>;

  // Constructor for a new, empty, cache element
  CacheEntry(const shape_t& default_shape) :
    chunk_data(default_shape), op_to_perform(CacheStatus::Nothing()) { }

  // Fill this cache element with data
  CacheEntry& operator=(std::tuple<const metadata_t&, const chunk_t&, const status_t&> other);
  
  metadata_t chunk_meta;
  chunk_t chunk_data;
  status_t op_to_perform;
};

// acts like a cached version of `NDVecArrayStreamer`
template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
class ChunkCache {

public:
  
  using chunk_t = ArrayT<T, dims, vec_dims>;
  using chunk_shape_t = typename chunk_t::shape_t;
  using chunk_meta_t = ChunkMetadata<dims>;

private:

  using cache_entry_t = CacheEntry<ArrayT, T, dims, vec_dims>;
  using status_t = typename cache_entry_t::status_t;
  
public:

  // `cache_size` ... number of chunks that can be kept in the cache
  // `init_cache_el_shape` ... initial shape to reserve for each element in the cache
  ChunkCache(std::filesystem::path workdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
	     const Vector<std::size_t, dims>& streamer_chunk_size,
	     std::size_t initial_buffer_size = 10000);

  // Makes some sensible simplified choices
  ChunkCache(std::filesystem::path workdir, std::size_t cache_size, std::size_t init_cache_el_linear_size = 400,
	     std::size_t initial_buffer_size = 10000);
  
  ~ChunkCache();

  // adds a new chunk with contents `chunk_data` and metadata `chunk_meta`
  // assumes that this is a new chunk and does not already exist
  void RegisterNewChunk(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data);

  // gets an existing chunk based on its metadata `chunk_meta`
  const chunk_t& RetrieveChunk(const chunk_meta_t& chunk_meta);

  // appends a slice `slice` to the existing chunk with metadata `chunk_meta`,
  // updates the metadata
  template <std::size_t axis>
  void AppendSlice(chunk_meta_t& chunk_meta, const chunk_t& slice);

  void FlushCache();
  
private:

  // deserialize chunk (and its metadata) from file and insert into the cache
  // returns the buffer slot of the inserted element
  cache_entry_t& deserialize_into_cache(const chunk_meta_t& chunk_meta);

  // inserts a new cache element into the cache
  void insert_into_cache(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data, const status_t& stat);
  
  // makes sure cache element with index `ind` is up-to date and ready for retrieval of values
  void sync_cache_element_for_read(cache_entry_t& cache_entry);

  void descope_cache_element(cache_entry_t& cache_entry);

  std::filesystem::path get_abs_path(const std::filesystem::path& chunk_path);
  
private:

  std::filesystem::path m_workdir;
  Vector<std::size_t, dims> m_streamer_chunk_size;
  Cache<std::size_t, cache_entry_t> m_cache;
  stor::NDVecArrayStreamer<ArrayT, T, dims, vec_dims> m_streamer;
  
};

// -------------------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
class ChunkLibrary {

public:

  using chunk_t = ArrayT<T, dims, vec_dims>;
  using chunk_shape_t = typename chunk_t::shape_t;
  using shape_t = typename ChunkIndex<dims>::shape_t;  
  using metadata_t = typename ChunkIndex<dims>::metadata_t;
  using view_t = typename chunk_t::view_t;  
  using ind_t = Vector<std::size_t, dims>;
  
public:

  ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size, std::size_t init_cache_el_linear_size = 400);
  ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
	       const Vector<std::size_t, dims>& streamer_chunk_size);

  // Add a new chunk
  void RegisterChunk(const ind_t& start_ind, const chunk_t& chunk);

  // Append a slice to an existing chunk along a certain axis
  template <std::size_t axis>
  void AppendSlice(const ind_t& start_ind, const chunk_t& slice);
  
  // Single-element retrieval
  view_t operator[](const ind_t& ind);

  // Other properties
  shape_t GetShape() {return m_index.GetShape();};

private:

  void CalculateShape();
  
private:

  std::filesystem::path m_libdir;
  std::filesystem::path m_index_path;

  ChunkIndex<dims> m_index;
  ChunkCache<ArrayT, T, dims, vec_dims> m_cache;
  
};

// ----------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
class DistributedNDVecArray {

public:

  using chunk_t = typename ChunkLibrary<ArrayT, T, dims, vec_dims>::chunk_t;
  using chunk_shape_t = typename ChunkLibrary<ArrayT, T, dims, vec_dims>::chunk_shape_t;
  using shape_t = typename ChunkLibrary<ArrayT, T, dims, vec_dims>::shape_t;
  using view_t = typename chunk_t::view_t;
  using ind_t = Vector<std::size_t, dims>;
  
public:

  DistributedNDVecArray(std::filesystem::path workdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
			const Vector<std::size_t, dims>& streamer_chunk_size);
  DistributedNDVecArray(std::filesystem::path workdir, std::size_t cache_size = 1, std::size_t init_cache_el_linear_size = 100);

  // Single-chunk operations
  void RegisterChunk(const ind_t& start_ind, const chunk_t& chunk);
  
  view_t operator[](const ind_t& ind);
  
  template <std::size_t axis>
  void AppendSlice(const ind_t& start_ind, const chunk_t& slice);

  shape_t GetShape(){ return m_library.GetShape(); };
  
  // Multi-chunk operations

  // Imports another distributed array and add it to this one
  void Import(const DistributedNDVecArray<ArrayT, T, dims, vec_dims>& other);
  
  // Rebalance chunks
  void RebuildChunks(const ind_t& requested_chunk_shape, std::filesystem::path tmpdir);
  
private:

  ChunkLibrary<ArrayT, T, dims, vec_dims> m_library;
  
};

#include "DistributedNDVecArray.hxx"
