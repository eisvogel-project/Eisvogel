#include <iostream>

#include "Cache.hh"
#include "NDVecArray.hh"
#include "DistributedNDVecArray.hh"

int main(int argc, char* argv[]) {

  auto print_element = [](VectorView<float, 2> elem) {
    std::cout << elem[0] << ", " << elem[1] << std::endl;
  };
  
  Vector<std::size_t, 3> shape{2u, 2u, 2u};
  using cache_entry_t = CacheEntry<NDVecArray, float, 3, 2>;


  float init_val = 0.0;
  Cache<std::size_t, cache_entry_t> testCache(4, shape, 0.0f);  

  std::cout << "inserting into cache" << std::endl;
  
  testCache.insert_no_overwrite(10, 123);
    
  // 

  // CacheEntry<NDVecArray, float, 3, 2> another_entry(shape, 10.0f);

  // std::cout << " ---- " << std::endl;

  // another_entry.chunk_data.loop_over_elements(print_element);
  
  // std::cout << " ---- " << std::endl;
  
  // another_entry = test_entry;

  // std::cout << " ---- " << std::endl;

  // another_entry.chunk_data.loop_over_elements(print_element);
  
  std::cout << "done" << std::endl;
}
