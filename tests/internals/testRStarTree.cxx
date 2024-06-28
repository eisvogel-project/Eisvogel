#include <iostream>
#include <format>
#include <random>

#include "RStarTree.hh"
#include "Vector.hh"
#include "IteratorUtils.hh"

using CoordT = std::size_t;
using PayloadT = float;
constexpr int dims = 3;
using TreeT = RStarTree<CoordT, dims, PayloadT>;

template <std::size_t dims>
float ind_sum(const Vector<std::size_t, dims>& ind) {

  std::mt19937 gen64;
  
  std::size_t hash = 0;
  std::size_t mulfact = 1;
  for(std::size_t i = 0; i < dims; i++) {
    hash += ind[i] * mulfact;
    mulfact *= 100;
  }
  gen64.seed(hash);
  
  std::size_t sum = 0;
  for(std::size_t i = 0; i < dims; i++) {
    sum += ind[i] * i;
  }

  float randval = static_cast <float> (gen64()) / static_cast <float> (gen64.max());
  return std::sin(sum) + randval;
}

template <std::size_t dims, class CallableT>
void fill_tree(TreeT& to_fill, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
	       const Vector<std::size_t, dims>& chunk_shape, CallableT&& filler) {

  auto tree_filler = [&](const Vector<std::size_t, dims>& chunk_start_ind, const Vector<std::size_t, dims>& chunk_end_ind) -> void {    
    PayloadT val = filler(chunk_start_ind);
    to_fill.InsertElement(val, chunk_start_ind, chunk_end_ind);
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, tree_filler);
}

template <std::size_t dims, class CallableT>
void check_tree(TreeT& to_check, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
		const Vector<std::size_t, dims>& chunk_shape, CallableT&& truth_value_calc) {

  std::size_t entries_checked = 0;
  
  auto tree_checker = [&](const Vector<std::size_t, dims>& chunk_start_ind, const Vector<std::size_t, dims>& chunk_end_ind) -> void {
    PayloadT truth_val = truth_value_calc(chunk_start_ind);
    PayloadT retrieved = to_check.Search(chunk_start_ind);
    assert(retrieved == truth_val);

    entries_checked++;
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, tree_checker);

  std::cout << "Visited " << entries_checked << " tree entries" << std::endl;
}

int main(void) {
    
  TreeT tree(100);

  Vector<CoordT, dims> start(0u);
  Vector<CoordT, dims> end(100u);
  Vector<CoordT, dims> chunk_shape(2u);

  auto tree_payload_calculator = [](const Vector<CoordT, dims>& coords) -> PayloadT {
    return ind_sum(coords);
  };

  // Put some entries into the tree
  fill_tree(tree, start, end, chunk_shape, tree_payload_calculator);

  // Retrieve the values and check them
  check_tree(tree, start, end, chunk_shape, tree_payload_calculator);
  
  // Vector<CoordT, dims> canvas_start{0u, 0u};
  // Vector<CoordT, dims> canvas_end{10000u, 10000u};
  // Vector<CoordT, dims> chunk_size{10u, 10u};
  
  // auto tree_adder = [&tree](const Vector<CoordT, dims>& chunk_start, const Vector<CoordT, dims>& chunk_end) -> void {
  //   Vector<CoordT, dims> margin(1u);    
  //   tree.InsertElement(3.14, chunk_start + margin, chunk_end - margin);
  // };  
  // IteratorUtils::index_loop_over_chunks(canvas_start, canvas_end, chunk_size, tree_adder);

  // tree.RebuildSTR();

  // tree.DumpJSONTreeStructure("./testtree_structure_final.yaml");

  // std::cout << "finished dumping" << std::endl;
  
  // // point query
  // PayloadT retrieved = tree.Search({5u, 5u});

  // std::cout << "retrieved: " << retrieved << std::endl;

  // // rectangle query
  // std::vector<std::reference_wrapper<const PayloadT>> res = tree.Search({0u, 0u}, {5u, 5u});

  // std::cout << "rectangle query gave " << res.size() << " results" << std::endl;
  // for(auto& cur_res: res) {
  //   std::cout << "retrieved: " << cur_res << std::endl;
  // }

  // tree.DumpJSONTreeStructure("./testtree_structure_final.yaml");
  
  std::cout << "done" << std::endl;
}
