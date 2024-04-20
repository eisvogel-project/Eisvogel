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

  // fully-specified constructor
  ChunkMetadata(const ChunkType& chunk_type, const std::filesystem::path& filepath, const id_t& chunk_id,
		const Vector<std::size_t, dims> start_ind, const Vector<std::size_t, dims> shape, std::size_t overlap);

  // pretty printing
  friend std::ostream& operator<< <dims> (std::ostream& stream, const ChunkMetadata& meta);
  
  // Accessor that grows the recorded shape of a chunk
  template <std::size_t axis>
  void GrowChunk(std::size_t shape_growth);

  // Accessor to swap two axes
  template <std::size_t axis_1, std::size_t axis_2>
  void SwapAxes();
  
  // data members
  
  // type of this chunk
  ChunkType chunk_type;
  
  // path to the file where this chunk lives
  std::filesystem::path filepath;

  // unique id: does not change when slices are appended
  id_t chunk_id;
  
  // offset within the file where this chunk lives
  std::streampos start_pos;

  // the shape of this chunk ...
  Vector<std::size_t, dims> shape;
  
  // ... as well as the indices of the elements marking the start and end element in the chunk
  Vector<std::size_t, dims> start_ind;
  Vector<std::size_t, dims> end_ind;

  // The number of samples (in each coordinate direction) this chunk overlaps with the neighbouring ones
  std::size_t overlap;

  // subtract this from a global index to get the corresponding index within the chunk array data structure
  Vector<std::size_t, dims> loc_ind_offset;
};

// -------------------------

// provides fast coordinate-indexed lookup of chunks
template <std::size_t dims>
class ChunkIndex {

public:

  using metadata_t = ChunkMetadata<dims>;
  using shape_t = Vector<std::size_t, dims>;
  using ind_t = Vector<std::size_t, dims>;
  
public:

  // build from stored index file
  ChunkIndex(std::filesystem::path index_path);
  ~ChunkIndex();

  metadata_t& RegisterChunk(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& shape, std::size_t overlap);
  
  // get chunk that contains the element with index `ind`
  const metadata_t& GetChunk(const Vector<std::size_t, dims>& ind) const;
  metadata_t& GetChunk(const Vector<std::size_t, dims>& ind);

  // get chunks that, taken together, cover the rectangular region between `start_ind` and `end_ind`
  std::vector<std::reference_wrapper<const metadata_t>> GetChunks(const Vector<std::size_t, dims>& start_ind,
								  const Vector<std::size_t, dims>& end_ind);  
  
  shape_t GetShape();

  // housekeeping and moving operations
  void FlushIndex();
  void MoveIndex(std::filesystem::path new_index_path);
  void ClearIndex();
  void ImportIndex(std::filesystem::path index_path);
  
  // iterators over chunk metadata sets
  auto begin();
  auto end();

public:
  
  // to check whether a specific index, or index range, is contained in or overlaps with a chunk
  static bool is_in_chunk(const metadata_t& chunk, const Vector<std::size_t, dims>& ind);
  static bool is_in_region(const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& shape,
			   const Vector<std::size_t, dims>& ind);  
  static bool chunk_overlaps_region(const metadata_t& chunk, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind);

  // calculate the actual start / end index of the region where `chunk` overlaps with a specified index range
  static Vector<std::size_t, dims> get_overlap_start_ind(const metadata_t& chunk, const Vector<std::size_t, dims>& start_ind);
  static Vector<std::size_t, dims> get_overlap_end_ind(const metadata_t& chunk, const Vector<std::size_t, dims>& end_ind);
  
public:

  Vector<std::size_t, dims> get_start_ind();
  Vector<std::size_t, dims> get_end_ind();
  
private:

  void invalidate_cached_index_metadata();
  void calculate_and_cache_index_metadata();
  Vector<std::size_t, dims> calculate_start_ind();
  Vector<std::size_t, dims> calculate_end_ind();
  
  metadata_t& find_chunk_by_index(const Vector<std::size_t, dims>& ind);   
  id_t get_next_chunk_id();
  
  void load_and_rebuild_index();
  std::vector<metadata_t> load_index_entries(std::filesystem::path index_path);
  
private:

  id_t m_next_chunk_id;
  std::filesystem::path m_index_path;

  // Overall metadata describing chunk index
  bool m_index_metadata_valid;
  shape_t m_shape;
  ind_t m_start_ind;
  ind_t m_end_ind;
  
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

  // Removes a chunk with given metadata `chunk_meta`; implies that this chunk cannot be retrieved again
  void RemoveChunk(const chunk_meta_t& chunk_meta);
  
  // Replace chunk with `chunk_meta` by a new chunk with `new_chunk_meta` as metadata with data at `new_chunk_data`
  void ReplaceChunk(const chunk_meta_t& chunk_meta, const chunk_meta_t& new_chunk_meta, const chunk_t& new_chunk_data);
  
  // gets an existing chunk based on its metadata `chunk_meta`
  const chunk_t& RetrieveChunk(const chunk_meta_t& chunk_meta);

  // appends a slice `slice` to the existing chunk with metadata `chunk_meta`,
  // updates the metadata
  template <std::size_t axis>
  void AppendSlice(chunk_meta_t& chunk_meta, const chunk_t& slice);

  // housekeeping operations
  void FlushCache();
  void MoveCache(std::filesystem::path new_workdir);
  void ClearCache();
  void ImportCache(std::filesystem::path workdir);
  
private:

  // deserialize chunk (and its metadata) from file and insert into the cache
  // returns the buffer slot of the inserted element
  cache_entry_t& deserialize_into_cache(const chunk_meta_t& chunk_meta);

  // inserts a new cache element into the cache
  void insert_into_cache(const chunk_meta_t& chunk_meta, const chunk_t& chunk_data, const status_t& stat);
  
  // makes sure cache element with index `ind` is synchronized with its disk image
  void sync_cache_element_with_disk(cache_entry_t& cache_entry);

  void descope_cache_element(cache_entry_t& cache_entry);

  std::filesystem::path get_abs_path(const std::filesystem::path& chunk_path);
  
private:

  std::filesystem::path m_workdir;
  Vector<std::size_t, dims> m_streamer_chunk_size;
  Cache<std::size_t, cache_entry_t> m_cache;
  std::size_t m_cache_size;
  stor::NDVecArrayStreamer<ArrayT, T, dims, vec_dims> m_streamer;

  static constexpr std::string_view m_suffix = ".chunk";
  
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

  ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size = 1, std::size_t init_cache_el_linear_size = 400);
  ChunkLibrary(std::filesystem::path libdir, std::size_t cache_size, const chunk_shape_t& init_cache_el_shape,
	       const Vector<std::size_t, dims>& streamer_chunk_size);

  // Add a new chunk to the libraryy
  void RegisterChunk(const chunk_t& chunk, const ind_t& start_ind, const ind_t& end_ind, std::size_t overlap);
  
  // Append a slice to an existing chunk along a certain axis
  template <std::size_t axis>
  void AppendSlice(const ind_t& start_ind, const chunk_t& slice);
  
  // Single-element retrieval
  view_t operator[](const ind_t& ind);
  
  // Retrieval of rectangular region
  void FillArray(chunk_t& array, const ind_t& input_start, const ind_t& input_end, const ind_t& output_start);

  // Change order of axes (along with corresponding change in memory layout)
  template <std::size_t axis_1, std::size_t axis_2>
  void SwapAxes();
  
  // Chunk-aware iterators over all elements
  template <class CallableT>
  constexpr void index_loop_over_elements(const ind_t& start_ind, const ind_t& end_ind, CallableT&& worker);

  template <class CallableT>
  constexpr void index_loop_over_elements(CallableT&& worker);
  
  // housekeeping
  void MoveLibrary(std::filesystem::path new_libdir);
  void ClearLibrary();
  void DeleteLibrary();
  void FlushLibrary();
  void ImportLibrary(std::filesystem::path libdir);
  
  // Other properties
  shape_t GetShape() {return m_index.GetShape();};
  ind_t GetStartInd() {return m_index.get_start_ind();};
  ind_t GetEndInd() {return m_index.get_end_ind();};
  std::filesystem::path GetLibdir() {return m_libdir;};
  
private:

  void CalculateShape();

  // Map a general `worker` over all chunks
  // The worker must have the signature (metadata_t& new_chunk_meta, const chunk_t& chunk_data, chunk_t& new_chunk_data)
  // and can modify the new metadata in-place (comes initialized to the current metadata)
  template <class CallableT>
  void map_over_chunks(CallableT&& worker);
  
public:

  static constexpr std::string_view m_index_filename = "index.bin";

private:
  
  std::filesystem::path m_libdir;
  std::filesystem::path m_index_path;

protected:
  
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
  void RegisterChunk(const chunk_t& chunk, const ind_t& start_ind);
  void RegisterChunk(const chunk_t& chunk, const ind_t& start_ind, const ind_t& end_ind, std::size_t overlap);
 
  template <std::size_t axis>
  void AppendSlice(const ind_t& start_ind, const chunk_t& slice);

  shape_t GetShape(){ return m_library.GetShape(); };

  // Single-element random access
  view_t operator[](const ind_t& ind);  

  // Moves this distributed array to a different location
  void Move(std::filesystem::path dest);
  
  // Imports another distributed array and add it to this one
  void Import(std::filesystem::path dir);
  
  // Rebalance chunks in-place, i.e. output will remain in the same location
  void RebuildChunks(const ind_t& requested_chunk_shape, std::filesystem::path tmpdir);
  
  template <class BoundaryCallableT>
  void RebuildChunks(const ind_t& requested_chunk_shape, std::filesystem::path tmpdir,
		     std::size_t overlap, BoundaryCallableT&& boundary_evaluator);

  // Rebalance chunks in a certain region and provide the output at `outdir`
  void RebuildChunksPartial(const ind_t& start_ind, const ind_t& end_ind,
			    const ind_t& requested_chunk_shape, std::filesystem::path outdir);

  // Rebuild chunks in a certain region while introducing overlap between chunks
  // Use `BoundaryCallableT` with signature (DistributedNDVecArray&, Vector<int, dims>&, )
  // to specify the value of `ind` that lies outside the shape of this distributed array
  template <class BoundaryCallableT>
  void RebuildChunksPartial(const ind_t& start_ind, const ind_t& end_ind,
			    const ind_t& requested_chunk_shape, std::filesystem::path outdir,
			    std::size_t overlap, BoundaryCallableT&& boundary_evaluator);
  
  // Change the order of the axes
  template <std::size_t axis_1, std::size_t axis_2>
  void SwapAxes();

  // Extract an array
  void FillArray(chunk_t& array, const ind_t& start_ind, const ind_t& end_ind, const ind_t& output_start);  
  
  // Chunk-aware iterators
  template <class CallableT>
  constexpr void index_loop_over_elements(CallableT&& worker);

  // Other housekeeping
  std::filesystem::path GetWorkdir() { return m_library.GetLibdir(); }
  void Flush() { m_library.FlushLibrary(); } 
  
private:

  // utility function to iterate over chunk indices of elements that sit outside the index range specified by `start_ind` and `end_ind`
  using signed_ind_t = Vector<int, dims>;
  template <class CallableT>
  static void index_loop_over_penetrating_chunk_elements(const ind_t& start_ind, const ind_t& end_ind,
							 const signed_ind_t& chunk_start_ind, const signed_ind_t& chunk_end_ind,
							 CallableT&& worker);
  
private:

  ChunkLibrary<ArrayT, T, dims, vec_dims> m_library;
  
};

// Type shortcuts for good semantics
using Distributed_RZT_ErEz_Array = DistributedNDVecArray<NDVecArray, scalar_t, 3, 2>;
using Distributed_TRZ_ErEz_Array = DistributedNDVecArray<NDVecArray, scalar_t, 3, 2>;

#include "DistributedNDVecArray.hxx"
