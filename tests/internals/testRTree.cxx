#include <iostream>

#include "RTree.hh"
#include "Vector.hh"

int main(void) {
  
  BoundingBox<std::size_t, 3> bbox_1({6u, 6u, 6u}, {10u, 10u, 10u});
  BoundingBox<std::size_t, 3> bbox_2({5u, 5u, 5u}, {11u, 12u, 13u});

  BoundingBox<std::size_t, 3> bbox_3(bbox_1, bbox_2);

  BoundingBox<std::size_t, 3> bbox_4({0u, 0u, 0u}, {6u, 7u, 10u});
  
  std::cout << bbox_3 << std::endl;
  std::cout << bbox_3.contains({5u, 5u, 4u}) << std::endl;
  std::cout << bbox_3.volume() << std::endl;

  std::cout << bbox_3.overlaps(bbox_4) << std::endl;
  std::cout << bbox_3.compute_overlap_volume(bbox_4) << std::endl;
  
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
  
  // constexpr int dims = 3;
  // using IndexT = Vector<std::size_t, dims>;
  // using PayloadT = float;

  // RTree<IndexT, PayloadT, dims> tree(100);

  // tree.InsertElement(3.14, {0u, 0u, 0u}, {10u, 10u, 10u});
  
  std::cout << "done" << std::endl;
}
