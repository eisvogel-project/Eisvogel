#include <iostream>
#include "Cache.hh"

int main(int argc, char* argv[]) {

  Cache<std::size_t, std::size_t> testCache(4);

  std::cout << "empty cache" << std::endl;
  testCache.print_old_to_new();

  testCache.insert_no_overwrite(13, 3);
  testCache.insert_no_overwrite(14, 4);
  testCache.insert_no_overwrite(15, 5);
  testCache.insert_no_overwrite(16, 6);
  
  std::cout << "after insert" << std::endl;
  testCache.print_old_to_new();

  std::size_t evicted_element = testCache.evict_oldest_full_cache();
  std::cout << "evicted element: " << evicted_element << std::endl;
  
  std::size_t index = 16;
  if(testCache.contains(index)) {
    std::cout << "element at " << index << ": " << testCache.get(index) << std::endl;

    std::size_t evicted_element = testCache.evict(index);
    std::cout << "evicted element: " << evicted_element << std::endl;

    std::cout << "after eviction" << std::endl;
    testCache.print_old_to_new();
  }
  
}
