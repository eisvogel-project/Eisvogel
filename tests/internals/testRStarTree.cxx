#include <iostream>
#include <format>
#include <random>
#include <chrono>

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

  std::cout << "Filling tree ... " << std::endl;

  std::size_t entries_added = 0;
  
  auto tree_filler = [&](const Vector<std::size_t, dims>& chunk_start_ind, const Vector<std::size_t, dims>& chunk_end_ind) -> void {    
    PayloadT val = filler(chunk_start_ind);
    PayloadT& inserted_val = to_fill.InsertElement(val, chunk_start_ind, chunk_end_ind);
    assert(inserted_val == val);
    
    entries_added++;
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, tree_filler);

  std::cout << "Added " << entries_added << " entries to tree" << std::endl;
}

template <std::size_t dims, class CallableT>
void check_tree_element_query(TreeT& to_check, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
			      const Vector<std::size_t, dims>& chunk_shape, CallableT&& truth_value_calc) {

  std::cout << "Testing tree element retrieval ..." << std::endl;
  
  std::size_t num_queries = 0;
  std::chrono::microseconds total_duration{0};
  
  auto element_checker = [&](const Vector<std::size_t, dims>& chunk_start_ind, const Vector<std::size_t, dims>& chunk_end_ind) -> void {
    PayloadT truth_val = truth_value_calc(chunk_start_ind);

    // Check that everything inside the chunk maps onto the correct value
    auto element_checker = [&](const Vector<std::size_t, dims>& cur_ind) -> void {      
      PayloadT retrieved = to_check.Search(cur_ind);
      assert(retrieved == truth_val);      
      num_queries++;
    };

    auto start = std::chrono::high_resolution_clock::now();
    IteratorUtils::index_loop_over_elements(chunk_start_ind, chunk_end_ind, element_checker);
    auto stop = std::chrono::high_resolution_clock::now();
    auto cur_duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);    
    total_duration += cur_duration;
  };
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, element_checker);

  std::cout << "Checked " << num_queries << " element queries in " << total_duration
	    << " -> " << (float)total_duration.count() / (float)num_queries << " us per query" << std::endl;
}

template <std::size_t dims, class CallableT>
void check_tree_range_query(TreeT& to_check, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
			    const Vector<std::size_t, dims>& chunk_shape, CallableT&& truth_value_calc) {

  std::cout << "Testing tree range retrieval ..." << std::endl;

  std::size_t num_queries = 0;
  
  auto range_checker = [&](const Vector<std::size_t, dims>& chunk_start_ind, const Vector<std::size_t, dims>& chunk_end_ind) -> void {
    PayloadT truth_val = truth_value_calc(chunk_start_ind);
    
    std::vector<std::reference_wrapper<const PayloadT>> retrieved = to_check.Search(chunk_start_ind, chunk_end_ind);

    // Expect to match with exactly one chunk
    assert(retrieved.size() == 1);    
    assert(retrieved[0] == truth_val);

    num_queries++;
  };

  auto start = std::chrono::high_resolution_clock::now();
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, range_checker);
  auto stop = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);    

  std::cout << "Checked " << num_queries << " range queries in " << total_duration
	    << " -> " << (float)total_duration.count() / (float)num_queries << " us per query" << std::endl;
}

int main(void) {
    
  TreeT tree(100);

  Vector<CoordT, dims> chunk_shape(2u);
  Vector<CoordT, dims> start(0u);
  Vector<CoordT, dims> end = chunk_shape * 100u;

  auto tree_payload_calculator = [](const Vector<CoordT, dims>& coords) -> PayloadT {
    return ind_sum(coords);
  };

  // Put some entries into the tree
  fill_tree(tree, start, end, chunk_shape, tree_payload_calculator);

  // Retrieve the values and check them
  check_tree_element_query(tree, start, end, chunk_shape, tree_payload_calculator);
  check_tree_range_query(tree, start, end, chunk_shape, tree_payload_calculator);

  // Rebuild the tree
  tree.RebuildSTR();

  // Check again
  check_tree_element_query(tree, start, end, chunk_shape, tree_payload_calculator);
  check_tree_range_query(tree, start, end, chunk_shape, tree_payload_calculator);
    
  std::cout << "done" << std::endl;
}
