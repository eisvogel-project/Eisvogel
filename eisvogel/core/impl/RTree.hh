#pragma once

#include <vector>
#include <cassert>
#include <numeric>

// Contiguous storage
template <class T>
class MemoryPool {

public:
  using IndT = std::size_t;
  
public:

  MemoryPool(std::size_t init_size);

  // Get reference to empty slot into which new element can be placed
  T& get_empty_slot();
  
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
template <class IndexT, typename PayloadIndT, typename NodeIndT, std::size_t MAX_NODESIZE>
struct Node : BoundingBox<IndexT> {
  
  // Default constructor
  Node();
  
  void mark_as_leaf(PayloadIndT payload_ind);
  void mark_as_internal();
  void add_child(NodeIndT child_ind);
  
  bool is_leaf;
  
  // If this is a leaf index, this is the pointer to the payload
  PayloadIndT payload_ind;
  
  // List of pointers to nodes that are children of this node
  std::size_t num_child_nodes;
  std::array<NodeIndT, MAX_NODESIZE> child_inds;
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

  using NodeIndT = MemoryPool<class NodeT>::IndT;
  using PayloadIndT = MemoryPool<PayloadT>::IndT;
  using TreeNode = Node<IndexT, PayloadIndT, NodeIndT, MAX_NODESIZE>;
  
  // Contiguous storage for all internal tree nodes (that define the structure of the tree)
  // and tree leaves (where the data lives)
  MemoryPool<TreeNode> m_nodes;
  MemoryPool<PayloadT> m_data;
};

#include "RTree.hxx"
