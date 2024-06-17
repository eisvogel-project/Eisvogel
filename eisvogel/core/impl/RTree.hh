#pragma once

#include <vector>
#include <cassert>
#include <numeric>

// Contiguous storage
template <class T, std::unsigned_integral IndT>
class MemoryPool {
  
public:

  MemoryPool(std::size_t init_size);

  // Get reference to empty slot into which new element can be placed, starting from the front or the back of the pool
  IndT get_empty_slot_front();
  IndT get_empty_slot_back();
  
  void free_slot(IndT slot_ind);
  
  // Element retrieval if the index of the slot is known
  T& operator[](IndT ind);
  const T& operator[](IndT ind) const;

  // Fills copies of the stored elements into the vector at `dest`
  void fill_elements(std::vector<T> dest);
  
  // Reset the memory pool
  void reset();

public:

  // debug output
  void print_free();
  void print_allocated();

  // Grow the memory pool to make room for new elements
  void grow();
  
private:

  void reset_slot_lists();
  
  // To test the status (allocated / free) of an element with a given index
  bool is_allocated(IndT slot_ind);
  bool is_free(IndT slot_ind);
  
private:

  struct Slot {

    static constexpr IndT INVALID = std::numeric_limits<IndT>::max();

    Slot() : data_ind(INVALID), next_data_ind(INVALID), prev_data_ind(INVALID) { };
    
    IndT data_ind;  // Index in the array `m_data` ("data index") where the data for this slot lives
    IndT next_data_ind;  // Points to the data index of the next slot of the same kind ("free" or "allocated")
    IndT prev_data_ind;  // Points to the data index of the previous slot of the same kind ("free" or "allocated")
  };
  
private:

  std::size_t m_init_size;

  // Use indices as relative pointers to access all elements; since the pool can grow, any absolute references
  // will get invalidated!
  IndT m_free_start;
  IndT m_free_end;

  IndT m_alloc_start;
  IndT m_alloc_end;
  
  std::vector<T> m_data;
  std::vector<Slot> m_slots;
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

// Tree entry
template <class IndexT, class PayloadT>
struct Entry : BoundingBox<IndexT> {
  PayloadT payload;
};

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE = 5, std::size_t MIN_NODESIZE = 2>
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
  using TreeEntry = Entry<IndexT, PayloadT>;

private:

  IndT choose_subtree(IndT start_node);
  
private:
  
  IndT m_root_node_ind;
  
  // Contiguous storage for all internal tree nodes (that define the structure of the tree)
  // and tree leaves (where the data lives)
  MemoryPool<TreeNode, IndT> m_nodes;
  MemoryPool<TreeEntry, IndT> m_entries;
};

#include "RTree.hxx"
