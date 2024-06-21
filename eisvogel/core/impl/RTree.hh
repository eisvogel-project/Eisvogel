#pragma once

#include <vector>
#include <cassert>
#include <numeric>
#include <fstream>
#include "Vector.hh"

// Contiguous storage
template <class T>
class MemoryPool {
  
public:

  MemoryPool(std::size_t init_size);

  // Get reference to empty slot into which new element can be placed, starting from the front of the pool
  std::size_t get_empty_slot_front(); 
  void free_slot(std::size_t slot_ind);
  
  // Element retrieval if the index of the slot is known
  T& operator[](std::size_t ind);
  const T& operator[](std::size_t ind) const;

  // Determine the slot where `obj` is seated
  std::size_t get_slot(const T& obj) const;
  
  // Fills copies of the stored elements into the vector at `dest`
  // void fill_elements(std::vector<T> dest);
  
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
  bool is_allocated(std::size_t slot_ind);
  bool is_free(std::size_t slot_ind);
  
private:

  struct Slot {

    static constexpr std::size_t INVALID = std::numeric_limits<std::size_t>::max();

    Slot() : next_data_ind(INVALID), prev_data_ind(INVALID) { };
    
    std::size_t next_data_ind;  // Points to the data index of the next slot of the same kind ("free" or "allocated")
    std::size_t prev_data_ind;  // Points to the data index of the previous slot of the same kind ("free" or "allocated")
  };

public:

  static constexpr std::size_t INVALID_SLOT = Slot::INVALID;
  
private:
  
  using SlotListDataT = std::vector<Slot>;

  bool is_valid_slot(std::size_t& ind);
  
  // Operations on doubly-linked lists:
  // A doubly-linked list is described by its `list_start` and `list_end`

  // Remove element with slot index `to_remove` from list
  void remove_from_slot_list(std::size_t& to_remove, std::size_t& list_start, std::size_t& list_end);

  // Add element with slot index `to_add` to the back of the list
  void add_to_slot_list_back(std::size_t& to_add, std::size_t& list_start, std::size_t& list_end);

  // Add element with slot index `to_add` to the beginning of the list
  void add_to_slot_list_front(std::size_t& to_add, std::size_t& list_start, std::size_t& list_end);

  // Helpers to iterate over linked list
  template <class CallableT>
  void loop_over_elements_front_to_back(std::size_t& list_start, CallableT&& worker);

  template <class CallableT>
  void loop_over_elements_back_to_front(std::size_t& list_end, CallableT&& worker);
  
private:

  std::size_t m_init_size;

  // Use indices as relative pointers to access all elements; since the pool can grow, any absolute references
  // will get invalidated!
  std::size_t m_free_start;
  std::size_t m_free_end;

  std::size_t m_alloc_start;
  std::size_t m_alloc_end;

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

  // The coordinates of the center point of this bounding box
  Vector<CoordT, dims> compute_center();
  
  // Distance between the center of this bounding box and the center of `bbox` ...
  CoordT compute_center_distance(const BoundingBox<CoordT, dims>& bbox);

  // ... or an explicitly-passed center point
  CoordT compute_center_distance(const Vector<CoordT, dims>& center);
  
  // Extends this bounding box (if needed) so that it also contains the passed bounding `bbox`
  void extend(const BoundingBox<CoordT, dims>& bbox);

  // Resets this bounding box
  void reset_bounding_box();

  // Pretty printing
  friend std::ostream& operator<< <CoordT, dims> (std::ostream& stream, const BoundingBox& bbox);
  
  // The start- and end coordinates describing this bounding box
  Vector<CoordT, dims> start_coords;
  Vector<CoordT, dims> end_coords;

  // Also keep track of the shape
  Vector<CoordT, dims> shape;
};

// =======================

// TODO: can think about moving MAX_NODESIZE and MIN_NODESIZE to constexpr float
template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE = 5, std::size_t MIN_NODESIZE = 2>
class RTree {

  static_assert(MAX_NODESIZE > 1);
  static_assert(MIN_NODESIZE < MAX_NODESIZE);

  static constexpr float REINSERT_P_FRAC = 0.30;

private:

  // Tree entry
  struct TreeEntry : BoundingBox<CoordT, dims> {
    
    // Construct an empty node with zero spatial extent and default-constructed payload
    TreeEntry();
    
    PayloadT payload;
  };

  // Tree node
  struct TreeNode : BoundingBox<CoordT, dims> {
    
  public:
    
    // Construct an empty node with zero spatial extent
    TreeNode();
    
    // To initialize as internal or leaf node
    void set_as_empty_leaf_node();
    void set_as_empty_internal_node();

    void add_child(std::size_t child_slot);
    
    bool is_leaf;
    std::size_t num_children;
    
    // List of pointers to nodes that are children of this node
    // These can either be other internal nodes (if `is_leaf == false`) or entries (if `is_leaf == true`)
    std::array<std::size_t, MAX_NODESIZE + 1> child_slots;
  };
  
public:

  // Constructs empty tree and reserves `init_slot_size` slots to hold elements
  RTree(std::size_t init_slot_size);

  // Insert a new element into the tree, given the element and the start- and end coordinates of its bounding box
  void InsertElement(const PayloadT& elem, const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords);

  // // Rebuild the tree and rebalance the nodes, if needed
  // void Rebuild();
  
  // const PayloadT& Search(const IndexT& ind) const;
  // PayloadT& Search(const IndexT& ind);
  // std::vector<std::reference_wrapper<const PayloadT&>> Search(const IndexT& start_ind, const IndexT& end_ind);
  
private:

  // Internal functions that do all the work
  // Note: nodes and entries are passed around through their "slots" in the respective memory pools. These slots are never
  // invalidated, even if the memory pool changes the actual memory allocation (which can happen during some of the recursive
  // tree rearrangements that are necessary).
  
  // Generators that build new leaves or entries
  std::size_t build_new_entry(const PayloadT& elem, const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords);
  std::size_t build_new_leaf_node(std::size_t entry_slot);
  std::size_t build_new_internal_node(const std::vector<std::size_t>& child_node_slots);

  // Insert the entry at `entry_slot` into the tree starting at the node at `start_node_slot`.
  // If any tree rearrangements occured during insertion, this will return the slot of the tree node that should be inserted
  // into the node at `node_slot`. This needs to be done by the caller.
  std::size_t insert(std::size_t entry_slot, std::size_t start_node_slot, bool first_insert = true);

  // Returns the slot of the child node of `start_node` where the entry at `entry_slot` should be inserted
  std::size_t choose_subtree(std::size_t entry_slot, std::size_t start_node_slot);
  
  // Handles an overfull node at `node_slot`; if as part of the clean-up an additional node was created, return its slot so that the caller can take care of it
  std::size_t overflow_treatment(std::size_t node_slot, bool first_insert);

  // Cleans up the tree starting at `node_slot` by reinserting its children
  void reinsert(std::size_t node_slot);

  // Cleans up the tree starting at `node_slot` by splitting the node and distributing its elements
  // Returns the slot of the newly-created node that needs to be added at the caller's level
  std::size_t split(std::size_t node_slot);
  
  // Fetch the bounding box of the child at `child_slot` of the node at `node_slot`
  BoundingBox<CoordT, dims>& get_bbox(std::size_t node_slot, std::size_t child_slot);

  // Force recalculation of the bounding box of the node at `node_slot`
  void recalculate_bbox(std::size_t node_slot);

  // Among all child nodes of the node at `node_slot`, find the one whose overlap with the other child nodes increases by the least amount
  // if the element with bounding box `to_add` was added to it, and return its node slot.
  // Resolves ties by returning the child that nees least volume enlargement.
  std::size_t find_min_overlap_enlargement_child(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add);

  // Among all child nodes of the node at `node_slot`, find the one whose volume increases by the least amount
  // if the element with bounding box `to_add` was added to it, and return its node slot.
  // Resolves ties by returning the child whose bounding box has the smallest volume.
  std::size_t find_min_volume_enlargement_child(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add);

  // Calculates the overlap enlargement (as defined by Beckmann et. al) of `child_slot` with all other children in the node at `node_slot`,
  // if `to_add` is added to the child at `child_slot`.
  std::size_t child_overlap_enlargement(std::size_t node_slot, std::size_t child_slot, const BoundingBox<CoordT, dims>& to_add);

  // Calculates the volume enlargement of the node at `node_slot` if `to_add` is added to `node_slot`
  std::size_t node_volume_enlargement(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add);

  // Gets the absolute volume of the node at `node_slot`
  std::size_t node_volume(std::size_t node_slot);
  
private:
  
  std::size_t m_root_slot;
  
  // Contiguous storage for all internal tree nodes (that define the structure of the tree)
  // and tree leaves (where the data lives)
  using NodePool = MemoryPool<TreeNode>;
  using EntryPool = MemoryPool<TreeEntry>;
  
  NodePool m_nodes;
  EntryPool m_entries;
};

#include "RTree.hxx"
