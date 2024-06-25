#include <iostream>

#include "RStarTree.hh"
#include "Vector.hh"
#include "IteratorUtils.hh"

int main(void) {
  
  using CoordT = std::size_t;
  constexpr int dims = 2;  // use a 2-dim tree so that visualizations become easier
  using PayloadT = float;
  
  RStarTree<CoordT, dims, PayloadT> tree(100);

  Vector<CoordT, dims> canvas_start{0u, 0u};
  Vector<CoordT, dims> canvas_end{100u, 100u};
  Vector<CoordT, dims> chunk_size{10u, 10u};

  auto tree_adder = [&tree](const Vector<CoordT, dims>& chunk_start, const Vector<CoordT, dims>& chunk_end) -> void {
    tree.InsertElement(3.14, chunk_start, chunk_end);
  };  
  IteratorUtils::index_loop_over_chunks(canvas_start, canvas_end, chunk_size, tree_adder);
  
  // point query
  PayloadT retrieved = tree.Search({5u, 5u});

  std::cout << "retrieved: " << retrieved << std::endl;

  // rectangle query
  std::vector<std::reference_wrapper<const PayloadT>> res = tree.Search({0u, 0u}, {5u, 5u});

  std::cout << "rectangle query gave " << res.size() << " results" << std::endl;
  for(auto& cur_res: res) {
    std::cout << "retrieved: " << cur_res << std::endl;
  }

  tree.DumpJSONTreeStructure("./testtree_structure.yaml");
  
  std::cout << "done" << std::endl;
}
