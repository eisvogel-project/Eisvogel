#include <iostream>
#include <fstream>
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

  ChunkIndex<3> chunk_index("index.bin");

  Vector<std::size_t, 3> start_ind{10u, 10u, 10u};
  Vector<std::size_t, 3> end_ind{13u, 13u, 13u};
  Vector<std::size_t, 3> chunk_shape = end_ind - start_ind;

  ChunkMetadata<3> chunk_meta(ChunkType::specified, "abc.bin", 3, start_ind, chunk_shape);

  std::fstream ofs;
  ofs.open("index.bin", std::ios::out | std::ios::binary);
  stor::Traits<ChunkMetadata<3>>::serialize(ofs, chunk_meta);
  ofs.close();

  std::fstream ifs;
  ifs.open("index.bin", std::ios::in | std::ios::binary);
  ChunkMetadata<3> chunk_meta_read = stor::Traits<ChunkMetadata<3>>::deserialize(ifs);
  ifs.close();

  std::cout << chunk_meta.filepath << std::endl;
  std::cout << chunk_meta_read.filepath << std::endl;
  
  std::cout << chunk_meta.start_ind[0] << ", " << chunk_meta.start_ind[1] << ", " << chunk_meta.start_ind[2] << std::endl;
  std::cout << chunk_meta_read.start_ind[0] << ", " << chunk_meta_read.start_ind[1] << ", " << chunk_meta_read.start_ind[2] << std::endl;
  
  std::cout << chunk_meta.end_ind[0] << ", " << chunk_meta.end_ind[1] << ", " << chunk_meta.end_ind[2] << std::endl;
  std::cout << chunk_meta_read.end_ind[0] << ", " << chunk_meta_read.end_ind[1] << ", " << chunk_meta_read.end_ind[2] << std::endl;
    
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
