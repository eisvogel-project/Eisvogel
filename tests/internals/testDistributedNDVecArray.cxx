#include <iostream>
#include <tuple>

#include "Cache.hh"
#include "NDVecArray.hh"
#include "DistributedNDVecArray.hh"

int main(int argc, char* argv[]) {

  auto print_element = [](VectorView<float, 2> elem) {
    std::cout << elem[0] << ", " << elem[1] << std::endl;
  };
  
  Vector<std::size_t, 3> shape{2u, 2u, 2u};
  Vector<std::size_t, 3> streamer_chunk_size{1u, stor::INFTY, stor::INFTY};
  
  ChunkCache<NDVecArray, float, 3, 2> chunk_cache(4, shape, streamer_chunk_size);

  ChunkMetadata<3> test_meta;

  chunk_cache.RetrieveChunk(test_meta);
  
  // using cache_entry_t = CacheEntry<NDVecArray, float, 3, 2>;


  
  // Cache<std::size_t, cache_entry_t> testCache(4, shape, 0.0f);  

  // std::cout << "creating local buffer" << std::endl;

  // ChunkMetadata<3> meta;
  // NDVecArray<float, 3, 2> buffer(shape, 1.0f);
  // CacheStatus stat = CacheStatus::nothing;
  
  // std::cout << "inserting into cache" << std::endl;

  // testCache.insert_no_overwrite(10, std::forward_as_tuple(meta, buffer, stat));
   

  // CacheEntry<NDVecArray, float, 3, 2> another_entry(shape, 10.0f);

  // std::cout << " ---- " << std::endl;

  // another_entry.chunk_data.loop_over_elements(print_element);
  
  // std::cout << " ---- " << std::endl;
  
  // another_entry = test_entry;

  // std::cout << " ---- " << std::endl;

  // another_entry.chunk_data.loop_over_elements(print_element);
  
  std::cout << "done" << std::endl;
}
