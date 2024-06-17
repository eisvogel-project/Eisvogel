#include <iostream>

#include "RTree.hh"
#include "Vector.hh"

int main(void) {

  MemoryPool<int> pool(10);
  
  constexpr int dims = 3;
  using IndexT = Vector<std::size_t, dims>;
  using PayloadT = float;

  RTree<IndexT, PayloadT, dims> tree(100);
  
  std::cout << "done" << std::endl;
}
