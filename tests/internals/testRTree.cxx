#include <iostream>

#include "RTree.hh"
#include "Vector.hh"

int main(void) {

  constexpr int dims = 3;
  using IndexT = Vector<std::size_t, dims>;
  using PayloadT = float;

  RTree<IndexT, PayloadT, dims> tree;
  
  std::cout << "done" << std::endl;
}
