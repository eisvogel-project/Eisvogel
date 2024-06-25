#include <iostream>

#include "RStarTree.hh"
#include "Vector.hh"

int main(void) {
  
  using CoordT = std::size_t;
  constexpr int dims = 3;
  using PayloadT = float;
  
  RStarTree<CoordT, dims, PayloadT> tree(100);
  tree.InsertElement(3.14, {0u, 0u, 0u}, {10u, 10u, 10u});

  // point query
  PayloadT retrieved = tree.Search({5u, 5u, 5u});

  std::cout << "retrieved: " << retrieved << std::endl;

  // rectangle query
  std::vector<std::reference_wrapper<const PayloadT>> res = tree.Search({0u, 0u, 0u}, {5u, 5u, 5u});

  std::cout << "rectangle query gave " << res.size() << " results" << std::endl;
  for(auto& cur_res: res) {
    std::cout << "retrieved: " << cur_res << std::endl;
  }

  tree.DumpJSONTreeStructure("./testtree_structure.yaml");
  
  std::cout << "done" << std::endl;
}
