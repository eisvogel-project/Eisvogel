#pragma once

#include <vector>
#include <cassert>
#include <numeric>

// Contiguous storage
template <class T, std::unsigned_integral IndT>
class MemoryPool {
  
public:

  MemoryPool(std::size_t init_size);

  // Get reference to empty slot into which new element can be placed
  IndT get_empty_slot();
  
  // Element retrieval
  T& operator[](IndT ind);
  const T& operator[](IndT ind) const;  

private:

  // Grow the memory pool to make room for new elements
  void grow();

  // To test the status (allocated / free) of an element with a given index
  bool is_allocated(IndT ind);
  bool is_free(IndT ind);
  
private:
  
  std::vector<T> m_data;
  std::vector<IndT> m_slots_free;
};

// =======================

// A general bounding box with start and end coordinates
template <class IndexT>
struct BoundingBox {
   
  IndexT start_ind;
  IndexT end_ind;
};

// Tree node
template <class IndexT, typename IndT, std::size_t MAX_NODESIZE>
struct Node : BoundingBox<IndexT> {
  
  // Default constructor
  Node();
  
  void mark_as_empty_leaf();
  void mark_as_empty_internal();
  void add_child(IndT child_ind);
  
  bool is_leaf;
  
  // List of pointers to nodes that are children of this node
  // These can either be other internal nodes (if `is_leaf == false`) or entries (if `is_leaf == true`)
  std::size_t num_child_nodes;
  std::array<IndT, MAX_NODESIZE> child_inds;
};

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE = 5>
class RTree {

public:

  // Constructs empty tree and reserves `init_slot_size` slots to hold elements
  RTree(std::size_t init_slot_size);

  void AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind);

  // Rebuild the tree and rebalance the nodes, if needed
  void Rebuild();
  
  const PayloadT& Search(const IndexT& ind) const;
  PayloadT& Search(const IndexT& ind);
  std::vector<std::reference_wrapper<const PayloadT&>> Search(const IndexT& start_ind, const IndexT& end_ind);

  
private:
  
  using IndT = std::size_t;
  using TreeNode = Node<IndexT, IndT, MAX_NODESIZE>;

  IndT m_root_node_ind;
  
  // Contiguous storage for all internal tree nodes (that define the structure of the tree)
  // and tree leaves (where the data lives)
  MemoryPool<TreeNode, IndT> m_nodes;
  MemoryPool<PayloadT, IndT> m_data;
};

#include "RTree.hxx"
