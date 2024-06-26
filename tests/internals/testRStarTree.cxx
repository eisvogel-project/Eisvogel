#include <iostream>

#include <format>

#include "RStarTree.hh"
#include "Vector.hh"
#include "IteratorUtils.hh"

int main(void) {
  
  using CoordT = std::size_t;
  constexpr int dims = 2;  // use a 2-dim tree so that visualizations become easier
  using PayloadT = float;
  
  RStarTree<CoordT, dims, PayloadT> tree(100);

  Vector<CoordT, dims> canvas_start{0u, 0u};
  Vector<CoordT, dims> canvas_end{10000u, 10000u};
  Vector<CoordT, dims> chunk_size{10u, 10u};
  
  auto tree_adder = [&tree](const Vector<CoordT, dims>& chunk_start, const Vector<CoordT, dims>& chunk_end) -> void {
    Vector<CoordT, dims> margin(1u);    
    tree.InsertElement(3.14, chunk_start + margin, chunk_end - margin);
  };  
  IteratorUtils::index_loop_over_chunks(canvas_start, canvas_end, chunk_size, tree_adder);

  tree.RebuildSTR();

  tree.DumpJSONTreeStructure("./testtree_structure_final.yaml");

  std::cout << "finished dumping" << std::endl;
  
  // point query
  PayloadT retrieved = tree.Search({5u, 5u});

  std::cout << "retrieved: " << retrieved << std::endl;

  // rectangle query
  std::vector<std::reference_wrapper<const PayloadT>> res = tree.Search({0u, 0u}, {5u, 5u});

  std::cout << "rectangle query gave " << res.size() << " results" << std::endl;
  for(auto& cur_res: res) {
    std::cout << "retrieved: " << cur_res << std::endl;
  }

  tree.DumpJSONTreeStructure("./testtree_structure_final.yaml");
  
  std::cout << "done" << std::endl;
}
