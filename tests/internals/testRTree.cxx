#include <iostream>

#include "RTree.hh"
#include "Vector.hh"

int main(void) {

  // using IndT = std::size_t;
  
  // MemoryPool<int, IndT> pool(3);
  // pool.print_free();
  // pool.print_allocated();

  // std::vector<IndT> allocated_slots;
  
  // for(std::size_t i = 0; i < 7; i++) {
  //   IndT slot = pool.get_empty_slot_front();
  //   allocated_slots.push_back(slot);
    
  //   std::cout << "allocated slot " << slot <<  std::endl;
  //   pool.print_free();
  //   pool.print_allocated();
  // }

  // for(std::size_t slot : allocated_slots) {
  //   std::cout << "freeing slot " << slot << std::endl;
  //   pool.free_slot(slot);
  //   pool.print_free();
  //   pool.print_allocated();
  // }
    
  // pool.grow();

  // pool.print_free();
  
  // pool[pool.get_empty_slot()] = 18;
  // std::cout << pool[9] << std::endl;
  
  constexpr int dims = 3;
  using IndexT = Vector<std::size_t, dims>;
  using PayloadT = float;

  RTree<IndexT, PayloadT, dims> tree(100);

  tree.InsertElement(3.14, {0u, 0u, 0u}, {10u, 10u, 10u});
  
  std::cout << "done" << std::endl;
}
