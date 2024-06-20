#pragma once

#include <vector>
#include <cassert>
#include <numeric>
#include <fstream>
#include "Vector.hh"

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

template <typename CoordT, std::size_t dims>
struct BoundingBox;

template <typename CoordT, std::size_t dims>
std::ostream& operator<<(std::ostream& stream, const BoundingBox<CoordT, dims>& bbox);

// A general bounding box with start and end coordinates
template <typename CoordT, std::size_t dims>
struct BoundingBox {

  // Default constructor
  BoundingBox();
  
  // Construct a new bounding box from start and end coordinates
  BoundingBox(const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords);

  // Copy constructor
  BoundingBox(const BoundingBox<CoordT, dims>& bbox);
  
  // Construct a new bounding box as the convex hull of two other bounding boxes `bbox_1` and `bbox_2`
  BoundingBox(const BoundingBox<CoordT, dims>& bbox_1, const BoundingBox<CoordT, dims>& bbox_2);
  
  // Checks if this bounding box contains the point with `coords`
  bool contains(const Vector<CoordT, dims>& coords);
  bool contains(const Vector<CoordT, dims>& coords) requires(std::same_as<CoordT, std::size_t>);

  // Checks if this bounding box overlaps with the passed `bbox`
  bool overlaps(const BoundingBox<CoordT, dims>& bbox);

  // Volume of this bounding box
  CoordT volume();

  // Overlap between this bounding box and `bbox`
  CoordT compute_overlap_volume(const BoundingBox<CoordT, dims>& bbox);
  
  // Extends this bounding box (if needed) so that it also contains the passed bounding `bbox`
  void extend(const BoundingBox<CoordT, dims>& bbox);

  // Pretty printing
  friend std::ostream& operator<< <CoordT, dims> (std::ostream& stream, const BoundingBox& bbox);
  
  // The start- and end coordinates describing this bounding box
  Vector<CoordT, dims> start_coords;
  Vector<CoordT, dims> end_coords;

  // Also keep track of the shape
  Vector<CoordT, dims> shape;
};

// =======================

// Tree node
template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
struct Node : BoundingBox<IndexT, dims> {

public:
  
  // Default constructor
  Node();
  
  void mark_as_empty_leaf();
  void mark_as_empty_internal();
  void add_child(SlotIndT child_ind);

  SlotIndT get_min_area_enlargement_child(SlotIndT& new_entry);
  SlotIndT get_min_overlap_enlargement_child(SlotIndT& new_entry);
  
  bool is_leaf;
  std::size_t num_child_nodes;

  // List of pointers to nodes that are children of this node
  // These can either be other internal nodes (if `is_leaf == false`) or entries (if `is_leaf == true`)
  std::array<SlotIndT, MAX_NODESIZE + 1> child_inds;
  
private:

  // Calculates the overlap of the child with index `child_ind` with all remaining
  // children of this node
  std::size_t get_child_overlap(SlotIndT& child_ind);
  
  BoundingBox<IndexT, dims>& get_child_bbox(SlotIndT& child_ind);
  
  // Returns references to the children in this node
  // std::vector<std::reference_wrapper<BoundingBox<IndexT, dims>>> dereference_children();  
  
  // Finds the index of the "minimum-element" child, where `comp` implements pairwise
  // comparison between two children
  template <class ComparatorT>
  SlotIndT min_child(ComparatorT&& comp);   
};

// Tree entry
template <class IndexT, std::size_t dims, class PayloadT>
struct Entry : BoundingBox<IndexT, dims> {

  // Default constructor
  Entry();
  
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

  // Internal insertion routine that is called recursively
  SlotIndT insert(const SlotIndT& entry_ind, const SlotIndT& node_ind, bool first_insert = true);

  // Returns the index of the tree node (leaf or internal) into which the entry with `entry_ind` should best be inserted
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
