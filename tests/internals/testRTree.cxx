#include <iostream>

#include "RTree.hh"
#include "Vector.hh"

int main(void) {

  MemoryPool<int, std::size_t> pool(10);
  pool.print_free();  
  
  pool.grow();

  pool.print_free();
  
  // pool[pool.get_empty_slot()] = 18;
  // std::cout << pool[9] << std::endl;
  
  // constexpr int dims = 3;
  // using IndexT = Vector<std::size_t, dims>;
  // using PayloadT = float;

  // RTree<IndexT, PayloadT, dims> tree(100);
  
  std::cout << "done" << std::endl;
}
