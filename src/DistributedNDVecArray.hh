#pragma once

#include <filesystem>
#include <fstream>
#include <tuple>
#include <vector>

#include "Cache.hh"
#include "Vector.hh"
#include "NDVecArrayStreamer.hh"

enum class ChunkType : std::size_t {
  specified = 0,
  all_null = 1
};

template <std::size_t dims>
struct ChunkMetadata {

  using id_t = std::size_t;

  // default constructor
  ChunkMetadata() :
    chunk_type(ChunkType::all_null), filepath(""), chunk_id(0), start_pos(0), start_ind(0), end_ind(0) { }
  
  // type of this chunk
  ChunkType chunk_type;
  
  // path to the file where this chunk lives
  std::filesystem::path filepath;

  // unique id: does not change when slices are appended
  id_t chunk_id;
  
  // offset within the file where this chunk lives
  std::streampos start_pos;

  // indices of the elements marking the start and end element in the chunk
  Vector<std::size_t, dims> start_ind;
  Vector<std::size_t, dims> end_ind;  
};

// provides fast coordinate-indexed lookup of chunks
template <template<std::size_t> class MetadataT,
	  std::size_t dims>
class ChunkIndex {
  
public:

  // build from stored index file
  ChunkIndex(std::filesystem::path index_path);

  // get chunk that contains the element with index `ind`
  ChunkMetadata<dims> GetChunk(const Vector<std::size_t, dims>& ind);

  // get chunks that, taken together, cover the rectangular region between `start_ind` and `end_ind`
  std::vector<ChunkMetadata<dims>> GetChunks(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind);
  
private:

  // for now, do simple linear search and remember to check last-accessed chunk first; use an R-tree later
    
};

// describes the status of this cache entry relative to the status on disk
// 2) it is encountered in a cache lookup: determines what needs to be done before the cache entry becomes valid for
//    the lookup
//    *) nothing: cache entry is up to date and elements can directly be retrieved
//    *) serialize: cache entry is up to date and elements can directly be retrieved
//    *) append: need to append the data and then re-read the full chunk before elements can be retrieved
enum class CacheStatus : std::size_t {
  nothing = 0,
  serialize = 1,
  append = 2
};

// The elements stored in the cache
template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
struct CacheEntry {  

  using metadata_t = ChunkMetadata<dims>;
  using chunk_t = ArrayT<T, dims, vec_dims>;
  using shape_t = typename chunk_t::shape_t;

  // Constructor for a new, empty, cache element
  CacheEntry(const shape_t& default_shape, const T& default_value) :
    chunk_data(default_shape, default_value), op_to_perform(CacheStatus::nothing) { }

  // Fill this cache element with data
  CacheEntry& operator=(std::tuple<metadata_t&, chunk_t&, CacheStatus&> other);
  
  metadata_t chunk_meta;
  chunk_t chunk_data;
  CacheStatus op_to_perform;

  void descope_entry();
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
  
public:

  // `cache_size` ... number of chunks that can be kept in the cache
  // `init_cache_el_shape` ... initial shape to reserve for each element in the cache
  ChunkCache(std::size_t cache_size, const chunk_shape_t& init_cache_el_shape);

  // adds a new chunk with contents `chunk_data` and metadata `chunk_meta`
  // assumes that this is a new chunk and does not already exist
  void RegisterNewChunk(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data);

  // gets an existing chunk based on its metadata `chunk_meta`
  const chunk_t& RetrieveChunk(const chunk_meta_t& chunk_meta);

  // appends a slice `slice` to the existing chunk with metadata `chunk_meta`,
  // returns metadata describing the new chunk resulting from the operation
  chunk_meta_t AppendSlice(const chunk_meta_t& chunk_meta, const chunk_t& slice);

private:

  // deserialize chunk (and its metadata) from file and insert into the cache
  // returns the buffer slot of the inserted element
  cache_entry_t& deserialize_into_cache(const chunk_meta_t& chunk_meta);

  // inserts a new cache element into the cache
  void insert_into_cache(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data, const CacheStatus& stat);
  
  // makes sure cache element with index `ind` is up-to date and ready for retrieval of values
  void sync_cache_element_for_read(cache_entry_t& cache_element);
  // {
  //    proceed according to `op_to_perform` found at this cache slot:
  //      0) `get` the element at the given cache slot (turns it into the most-recently accessed one)
  //      1) `nothing`: do nothing
  //      2) `serialize`: this is already the full chunk, do nothing
  //      3) `append`: perform appending to disk, then reread into this slot (using the reference obtained at step 0)
  // }
  
private:

  stor::NDVecArrayStreamer<ArrayT, T, dims, vec_dims> m_streamer;
  Cache<std::size_t, cache_entry_t> m_cache;
  
};

template <typename T, std::size_t dims, std::size_t vec_dims>
class ChunkLibrary {

public:

  // single-element retrieval:
  // -> check if element is in cache -> if so, return
  // -> if not in cache, load containing chunk from disk and then return

  // element range retrieval:
  // -> get list of chunks that have some intersection with the requested range from the full chunk index
  // -> load them one-by-one and take out the intersection of the loaded cache with the requested range
    
private:

  // load into cache:
  // -> get cache insertion position
  // -> sync contained element to disk if needed
  // -> load requested chunk into freed-up buffer location

  // ----
  
  // data members
  
  // streamer
  // chunk cache
  // chunk index
  
};

// template <typename T, std::size_t dims, std::size_t vec_dims>
// class DistributedNDVecArray {

//   // handles all the multi-chunk stuff by talking to ChunkLibrary
//   // temporarily, can also have more than one chunk library (e.g. when rechunking)
  
// public:

//   DistributedNDVecArray(std::filesystem::path dirpath, std::size_t cache_size);

// private:

//   const std::filesystem::path m_dirpath;

  
//   // chunk library
// };

#include "DistributedNDVecArray.hxx"
