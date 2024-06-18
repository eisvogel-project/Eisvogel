#pragma once

#include <vector>
#include <cassert>
#include <numeric>

// Contiguous storage
template <class T, std::unsigned_integral SlotIndT>
class MemoryPool {
  
public:

  MemoryPool(std::size_t init_size);

  // Get reference to empty slot into which new element can be placed, starting from the front of the pool
  SlotIndT get_empty_slot_front(); 
  void free_slot(SlotIndT slot_ind);
  
  // Element retrieval if the index of the slot is known
  T& operator[](SlotIndT ind);
  const T& operator[](SlotIndT ind) const;

  // Fills copies of the stored elements into the vector at `dest`
  void fill_elements(std::vector<T> dest);
  
  // Reset the memory pool
  void reset();
  
public:

  // debug output: mark as private after testing
  void print_free();
  void print_allocated();

  // Grow the memory pool to make room for new elements
  void grow(); // mark as private after testing
  
private:

  void reset_slot_lists();
  
  // To test the status (allocated / free) of an element with a given index
  bool is_allocated(SlotIndT slot_ind);
  bool is_free(SlotIndT slot_ind);
  
private:

  struct Slot {

    static constexpr SlotIndT INVALID = std::numeric_limits<SlotIndT>::max();

    Slot() : next_data_ind(INVALID), prev_data_ind(INVALID) { };
    
    SlotIndT next_data_ind;  // Points to the data index of the next slot of the same kind ("free" or "allocated")
    SlotIndT prev_data_ind;  // Points to the data index of the previous slot of the same kind ("free" or "allocated")
  };

public:

  static constexpr SlotIndT INVALID_SLOT = Slot::INVALID;
  
private:
  
  using SlotListDataT = std::vector<Slot>;

  bool is_valid_slot(SlotIndT& ind);
  
  // Operations on doubly-linked lists:
  // A doubly-linked list is described by its `list_start` and `list_end`

  // Remove element with slot index `to_remove` from list
  void remove_from_slot_list(SlotIndT& to_remove, SlotIndT& list_start, SlotIndT& list_end);

  // Add element with slot index `to_add` to the back of the list
  void add_to_slot_list_back(SlotIndT& to_add, SlotIndT& list_start, SlotIndT& list_end);

  // Add element with slot index `to_add` to the beginning of the list
  void add_to_slot_list_front(SlotIndT& to_add, SlotIndT& list_start, SlotIndT& list_end);

  // Helpers to iterate over linked list
  template <class CallableT>
  void loop_over_elements_front_to_back(SlotIndT& list_start, CallableT&& worker);

  template <class CallableT>
  void loop_over_elements_back_to_front(SlotIndT& list_end, CallableT&& worker);
  
private:

  std::size_t m_init_size;

  // Use indices as relative pointers to access all elements; since the pool can grow, any absolute references
  // will get invalidated!
  SlotIndT m_free_start;
  SlotIndT m_free_end;

  SlotIndT m_alloc_start;
  SlotIndT m_alloc_end;

  std::vector<T> m_data;      // Memory pool has contiguous storage
  SlotListDataT m_slots;  // Metadata to keep track of the individual slots
};

// =======================

// A general bounding box with start and end coordinates
template <class IndexT, std::size_t dims>  // TODO: put a `requires` constraint so that `IndexT` must implement the [] operator
struct BoundingBox {

  // Checks if this bounding box contains the point with `ind`
  bool contains(const IndexT& ind);

  // Checks if this bounding box overlaps with the passed `bbox`
  bool overlaps(const BoundingBox& bbox);

  // Stretches this bounding box (if needed) so that it also contains the passed bounding `bbox`
  void stretch(const BoundingBox& bbox);
  
  IndexT start_ind;
  IndexT end_ind;
};

// Tree node
template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
struct Node : BoundingBox<IndexT, dims> {
  
  // Default constructor
  Node();
  
  void mark_as_empty_leaf();
  void mark_as_empty_internal();
  void add_child(SlotIndT child_ind);  
  
  bool is_leaf;
  
  // List of pointers to nodes that are children of this node
  // These can either be other internal nodes (if `is_leaf == false`) or entries (if `is_leaf == true`)
  std::size_t num_child_nodes;
  std::array<SlotIndT, MAX_NODESIZE + 1> child_inds;
};

// Tree entry
template <class IndexT, std::size_t dims, class PayloadT>
struct Entry : BoundingBox<IndexT, dims> {
  PayloadT payload;
};

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE = 5, std::size_t MIN_NODESIZE = 2>
class RTree {

  static_assert(MAX_NODESIZE > 1);
  static_assert(MIN_NODESIZE < MAX_NODESIZE);
  
public:

  // Constructs empty tree and reserves `init_slot_size` slots to hold elements
  RTree(std::size_t init_slot_size);

  void InsertElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind);

  // Rebuild the tree and rebalance the nodes, if needed
  void Rebuild();
  
  const PayloadT& Search(const IndexT& ind) const;
  PayloadT& Search(const IndexT& ind);
  std::vector<std::reference_wrapper<const PayloadT&>> Search(const IndexT& start_ind, const IndexT& end_ind);
  
private:

  using SlotIndT = std::size_t;
  using TreeNode = Node<IndexT, dims, SlotIndT, MAX_NODESIZE>;
  using TreeEntry = Entry<IndexT, dims, PayloadT>;
  
  using NodePool = MemoryPool<TreeNode, SlotIndT>;
  using EntryPool = MemoryPool<TreeEntry, SlotIndT>;

private:

  SlotIndT insert(const SlotIndT& entry_ind, const SlotIndT& node_ind, bool first_insert = true);
  SlotIndT choose_subtree(const SlotIndT& start_node, const SlotIndT& entry_ind);
  SlotIndT overflow_treatment(const SlotIndT& node_ind, bool first_insert);
  
private:
  
  SlotIndT m_root_node_ind;
  
  // Contiguous storage for all internal tree nodes (that define the structure of the tree)
  // and tree leaves (where the data lives)
  NodePool m_nodes;
  EntryPool m_entries;
};

#include "RTree.hxx"
