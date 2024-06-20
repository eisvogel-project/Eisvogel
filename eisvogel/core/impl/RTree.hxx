#include <algorithm>

template <class T>
MemoryPool<T>::MemoryPool(std::size_t init_size) : m_init_size(init_size),
						   m_free_start(Slot::INVALID), m_free_end(Slot::INVALID),
						   m_alloc_start(Slot::INVALID), m_alloc_end(Slot::INVALID),
						   m_data(init_size), m_slots(init_size) { 
  reset_slot_lists();
}

template <class T>
void MemoryPool<T>::reset() {

  m_data.clear();
  m_data.resize(m_init_size);

  m_slots.clear();
  m_slots.resize(m_init_size);
  
  reset_slot_lists();
}

template <class T>
void MemoryPool<T>::reset_slot_lists() {

  assert(m_data.size() == m_slots.size());

  // Reset the two doubly-linked lists
  m_alloc_start = Slot::INVALID;
  m_alloc_end = Slot::INVALID;

  m_free_start = Slot::INVALID;
  m_free_end = Slot::INVALID;
  
  // Mark everything as free
  for(std::size_t i = 0; i < m_slots.size(); i++) {
    add_to_slot_list_back(i, m_free_start, m_free_end);
  }
}

template <class T>
T& MemoryPool<T>::operator[](std::size_t ind) {
  return m_data[ind];
}

template <class T>
const T& MemoryPool<T>::operator[](std::size_t ind) const {
  return m_data[ind];
}

template <class T>
std::size_t MemoryPool<T>::get_slot(const T& obj) const {
  ptrdiff_t diff = &obj - &m_data.front();
  assert((diff >= 0) && (diff < m_data.size()));
  return diff;
}

template <class T>
std::size_t MemoryPool<T>::get_empty_slot_front() {

  // Grow the pool if needed
  if(m_free_start == Slot::INVALID) {
    grow();
  }

  // Get the index of the next free slot
  std::size_t alloc_slot_ind = m_free_start;  
  
  // Remove the newly allocated slot from the beginning of the list of free slots ...
  remove_from_slot_list(alloc_slot_ind, m_free_start, m_free_end);

  // ... and add it to the back of the list of allocated slots
  add_to_slot_list_back(alloc_slot_ind, m_alloc_start, m_alloc_end);

  // Return the allocated slot
  return alloc_slot_ind;
}

template <class T>
void MemoryPool<T>::free_slot(std::size_t ind) {

  // Remove the slot with index `ind` from the list of allocated slots ...
  remove_from_slot_list(ind, m_alloc_start, m_alloc_end);

  // ... and add it to the back of the list of free slots
  add_to_slot_list_back(ind, m_free_start, m_free_end);
}

template <class T>
void MemoryPool<T>::grow() {

  // Double the size of the currently-allocated memory
  std::size_t current_size = m_data.size();
  std::size_t new_size = 2 * current_size;
  m_data.resize(new_size);

  // Extend the slot list
  std::vector<Slot> new_slots_free(new_size - current_size);
  m_slots.insert(m_slots.end(), new_slots_free.begin(), new_slots_free.end());

  // Mark the new slots as free
  for(std::size_t i = current_size; i < new_size; i++) {
    add_to_slot_list_back(i, m_free_start, m_free_end);
  }  
}

template <class T>
bool MemoryPool<T>::is_valid_slot(std::size_t& ind) {
  return (ind != Slot::INVALID) && (ind < m_slots.size());
}

template <class T>
void MemoryPool<T>::remove_from_slot_list(std::size_t& to_remove, std::size_t& list_start, std::size_t& list_end) {

  assert(is_valid_slot(to_remove));
  
  // Get the indices of the previous and next element from the list
  std::size_t& next_elem = m_slots[to_remove].next_data_ind;
  std::size_t& prev_elem = m_slots[to_remove].prev_data_ind;

  if(prev_elem != Slot::INVALID) {

    // Connect the "previous" element to the "next" element
    m_slots[prev_elem].next_data_ind = next_elem;    
  }
  else {

    // The removed element was the first entry in the list; update `list_start`
    list_start = next_elem;
  }

  if(next_elem != Slot::INVALID) {

    // Connect the "next" element to the "previous" one
    m_slots[next_elem].prev_data_ind = prev_elem;
  }
  else {

    // The removed element was the last entry in the list; update `list_end`
    list_end = prev_elem;
  }
}

template <class T>
void MemoryPool<T>::add_to_slot_list_back(std::size_t& to_add, std::size_t& list_start, std::size_t& list_end) {

  assert(is_valid_slot(to_add));

  if(list_end != Slot::INVALID) {

    // Non-empty list: link the last element to the newly-added one
    m_slots[list_end].next_data_ind = to_add;
    m_slots[to_add].prev_data_ind = list_end;
    m_slots[to_add].next_data_ind = Slot::INVALID;
    list_end = to_add;
  }
  else {

    // Consistency: for an empty list, both `list_start` and `list_end` must be invalid
    assert(list_start == Slot::INVALID);
    
    // Empty list: `to_add` is the first element in the list
    list_start = to_add;
    list_end = to_add;
    m_slots[to_add].prev_data_ind = Slot::INVALID;
    m_slots[to_add].next_data_ind = Slot::INVALID;
  }
}

template <class T>
void MemoryPool<T>::add_to_slot_list_front(std::size_t& to_add, std::size_t& list_start, std::size_t& list_end) {

  assert(is_valid_slot(to_add));

  if(list_start != Slot::INVALID) {

    // Non-empty list: list the first element to the newly-added one
    m_slots[list_start].prev_data_ind = to_add;
    m_slots[to_add].next_data_ind = list_start;
    m_slots[to_add].prev_data_ind = Slot::INVALID;
    list_start = to_add;
  }
  else {

    // Consistency: for an empty list, both `list_start` and `list_end` must be invalid
    assert(list_end == Slot::INVALID);

    // Empty list: `to_add` is the first element in the list
    list_start = to_add;
    list_end = to_add;
    m_slots[to_add].prev_data_ind = Slot::INVALID;
    m_slots[to_add].next_data_ind = Slot::INVALID;
  }
}

template <class T>
template <class CallableT>
void MemoryPool<T>::loop_over_elements_front_to_back(std::size_t& list_start, CallableT&& worker) {

  std::size_t cur_elem = list_start;
  while(cur_elem != Slot::INVALID) {
    
    // call worker
    worker(cur_elem);
    // worker(m_data[cur_elem]);

    // Walk to the following element in the list
    cur_elem = m_slots[cur_elem].next_data_ind;
  }  
}

template <class T>
template <class CallableT>
void MemoryPool<T>::loop_over_elements_back_to_front(std::size_t& list_end, CallableT&& worker) {

  std::size_t cur_elem = list_end;
  while(cur_elem != Slot::INVALID) {

    // call worker
    worker(cur_elem);
    // worker(m_data[cur_elem]);

    // Walk to the previoud element in the list
    cur_elem = m_slots[cur_elem].prev_data_ind;
  }
}

// ------

template <class T>
void MemoryPool<T>::print_free() {

  auto printer = [](const std::size_t& cur_slot) {
    std::cout << cur_slot << " ";
  };

  std::cout << "----------------- " << std::endl;
  std::cout << "Free:" << std::endl;
  std::cout << "----------------- " << std::endl;
  loop_over_elements_front_to_back(m_free_start, printer);
  std::cout << std::endl;

  loop_over_elements_back_to_front(m_free_end, printer);
  std::cout << std::endl;
  std::cout << "----------------- " << std::endl;
}

template <class T>
void MemoryPool<T>::print_allocated() {

  auto printer = [](const std::size_t& cur_slot) {
    std::cout << cur_slot << " ";
  };

  std::cout << "----------------- " << std::endl;
  std::cout << "Allocated:" << std::endl;
  loop_over_elements_front_to_back(m_alloc_start, printer);
  std::cout << std::endl;

  loop_over_elements_back_to_front(m_alloc_end, printer);
  std::cout << std::endl;
  std::cout << "----------------- " << std::endl;
}

// ------

template <class T>
bool MemoryPool<T>::is_free(std::size_t ind) {
  return false; // TODO
}

template <class T>
bool MemoryPool<T>::is_allocated(std::size_t ind) {
  return !is_free(ind);
}

// =======================

template <class CoordT, std::size_t dims>
BoundingBox<CoordT, dims>::BoundingBox(const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) :
  start_coords(start_coords), end_coords(end_coords), shape(end_coords - start_coords) { }

template <class CoordT, std::size_t dims>
BoundingBox<CoordT, dims>::BoundingBox(const BoundingBox<CoordT, dims>& bbox) : BoundingBox<CoordT, dims>(bbox.start_coords, bbox.end_coords) { }

template <class CoordT, std::size_t dims>
BoundingBox<CoordT, dims>::BoundingBox(const BoundingBox<CoordT, dims>& bbox_1, const BoundingBox<CoordT, dims>& bbox_2) : BoundingBox<CoordT, dims>(bbox_1) {
  extend(bbox_2);
}

template <class CoordT, std::size_t dims>
void BoundingBox<CoordT, dims>::extend(const BoundingBox<CoordT, dims>& bbox) {

  // To take the convex hull, take the elementwise minimum or -maximum of the two start- or end coordinates ...
  start_coords = VectorUtils::min(start_coords, bbox.start_coords);
  end_coords = VectorUtils::max(end_coords, bbox.end_coords);

  // ... and update the shape
  shape = end_coords - start_coords;
}

template <class CoordT, std::size_t dims>
bool BoundingBox<CoordT, dims>::contains(const Vector<CoordT, dims>& coords) {
  throw std::logic_error("Not implemented yet!");
}

template <class CoordT, std::size_t dims>
bool BoundingBox<CoordT, dims>::contains(const Vector<CoordT, dims>& coords) requires(std::same_as<CoordT, std::size_t>) {

  for(std::size_t i = 0; i < dims; i++) {
    
    // Efficient out-of-range comparison with a single branch
    if((coords[i] - start_coords[i]) >= shape[i]) {
      return false;
    }
  }  
  return true;
}

template <class CoordT, std::size_t dims>
bool BoundingBox<CoordT, dims>::overlaps(const BoundingBox<CoordT, dims>& bbox) {

  for(std::size_t i = 0; i < dims; i++) {

    // the `start_ind` of the specified region must lie "to the left of" the chunk endpoint in all directions
    if(start_coords[i] >= bbox.end_coords[i]) {
      return false;
    }

    // the `end_ind` of the specified region must lie "to the right of" the chunk endpoint in all directions
    if(end_coords[i] <= bbox.start_coords[i]) {
      return false;
    }    
  }
  
  return true;
}

template <class CoordT, std::size_t dims>
CoordT BoundingBox<CoordT, dims>::volume() {
  return std::accumulate(shape.cbegin(), shape.cend(), (CoordT)1u, std::multiplies<CoordT>());
}

template <class CoordT, std::size_t dims>
CoordT BoundingBox<CoordT, dims>::compute_overlap_volume(const BoundingBox<CoordT, dims>& bbox) {

  // Start- and end coordinates of the overlapping region
  Vector<CoordT, dims> overlap_start_coords = VectorUtils::max(start_coords, bbox.start_coords);
  Vector<CoordT, dims> overlap_end_coords = VectorUtils::min(end_coords, bbox.end_coords);
    
  CoordT overlap_volume = (CoordT)1u;
  for(std::size_t i = 0; i < dims; i++) {

    // If the overlap has non-negative extent along a particular direction ...
    if(overlap_end_coords[i] > overlap_start_coords[i]) {

      // ... take it into account in its volume
      overlap_volume *= (overlap_end_coords[i] - overlap_start_coords[i]);
    }
    else {

      // Otherwise there is no finite overlap
      return (CoordT)0u;
    }
  }  
  return overlap_volume;
}

template <class CoordT, std::size_t dims>
std::ostream& operator<<(std::ostream& stream, const BoundingBox<CoordT, dims>& bbox) {

  // TODO: use box-drawing characters to make this nicer
  stream << "-----------------------------------------\n";
  stream << "start_coords: " << bbox.start_coords << "\n";
  stream << "end_coords: " << bbox.end_coords     << "\n";
  stream << "-----------------------------------------\n";

  return stream;
}

// =======================

// Methods to compare different bounding boxes (needed to determine optimal insertion position)

// template <class BoundedT>
// struct CompareByVolumeEnlargement {

//   CompareByVolumeEnlargement(const BoundedT& new_element) : m_new_element(new_element) { };

//   // Compares `bd_1` and `bd_2`: returns `true` if adding `new_element` to `bd_1` results
//   // in _less_ volume enlargement than adding `new_element` to `bd_2`
//   bool operator()(const BoundedT& bd_1, const BoundedT& bd_2) {
//     return true;
//   }

// private:
//   const BoundedT& m_new_element;
// };

// template <class BoundedT>
// struct CompareByOverlapEnlargement {

//   CompareByOverlapEnlargement(const BoundedT& new_element, std::vector<std::reference_wrapper<BoundedT>>& other_elements) :
//     m_new_element(new_element), m_other_elements(other_elements) { };

//   // Compares `bd_1` and `bd_2`: returns `true` if adding `new_element` to `bd_1` results
//   // in less _overlap_ enlargement than adding `new_element` to `bd_2`
//   // (The overlap is calculated with respect to the `other_elements`.)
//   bool operator()(const BoundedT& bd_1, const BoundedT& bd_2) {
//     return true;
//   }
  
// private:
//   const BoundedT& m_new_element;
//   const std::vector<std::reference_wrapper<BoundedT>> m_other_elements;
// };

// =======================

// =======================





// template <typename CoordT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
// SlotIndT Node<CoordT, dims, SlotIndT, MAX_NODESIZE>::get_min_area_enlargement_child(SlotIndT& new_entry_ind) {

//   // Does not make sense to call this function if there are no child nodes
//   assert(num_child_nodes > 0);

//   SlotIndT min_child_ind;
//   for(std::size_t i = 0; i < num_child_nodes; i++) {
    
//   }
  
//   return min_child_ind;
// }

// template <typename CoordT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
// std::size_t Node<CoordT, dims, SlotIndT, MAX_NODESIZE>::get_child_overlap(SlotIndT& child_ind) {

//   return 0;
// }

// template <typename CoordT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
// SlotIndT Node<CoordT, dims, SlotIndT, MAX_NODESIZE>::get_min_overlap_enlargement_child(SlotIndT& new_entry) {
  
//   return 0;
// }

// template <typename CoordT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
// template <class ComparatorT>
// SlotIndT Node<CoordT, dims, SlotIndT, MAX_NODESIZE>::min_child(ComparatorT&& comp) {
  
//   return 0;
// }

// ===============

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeEntry::TreeEntry() : BoundingBox<CoordT, dims>(0, 0), payload() { }

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeNode::TreeNode() : BoundingBox<CoordT, dims>(0, 0), is_leaf(true), num_children(0) {
  child_slots.fill(0);
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeNode::set_as_empty_leaf_node() {
  is_leaf = true;
  num_children = 0;
  child_slots.fill(0);  
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeNode::set_as_empty_internal_node() {
  is_leaf = false;
  num_children = 0;
  child_slots.fill(0);
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeNode::add_child(std::size_t child_slot) {
  assert(num_children <= MAX_NODESIZE);
  child_slots[num_children] = child_slot;
  num_children++;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::RTree(std::size_t init_slot_size) : m_root_slot(MemoryPool<TreeNode>::INVALID_SLOT),
											       m_nodes(init_slot_size), m_entries(init_slot_size) { }

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::build_new_entry(const PayloadT& elem,
										       const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) {
  std::size_t entry_slot = m_entries.get_empty_slot_front();
  TreeEntry& entry = m_entries[entry_slot];
  
  entry.payload = elem;
  entry.start_coords = start_coords;
  entry.end_coords = end_coords;
  
  return entry_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::build_new_leaf_node(std::size_t entry_slot) {
  
  std::size_t node_slot = m_nodes.get_empty_slot_front();
  TreeNode& node = m_nodes[node_slot];

  // The new node is a leaf ...
  node.set_as_empty_leaf_node();
  
  // ... and holds a reference to the new entry
  node.add_child(entry_slot);
  
  return node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::InsertElement(const PayloadT& elem,
									      const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) {

  // Take ownership of the `elem` and store it as an entry in the memory pool
  std::size_t entry_slot = build_new_entry(elem, start_coords, end_coords);

  std::cout << "put new element into slot " << entry_slot << std::endl;

  // Insert the entry into the node structure of the tree
  if(m_root_slot == NodePool::INVALID_SLOT) {

    // Empty tree: build new root node
    m_root_slot = build_new_leaf_node(entry_slot);
  }
  else {

    // Non-empty tree: traverse the tree starting from the root node and insert the new entry at the correct location
    insert(entry_slot, m_root_slot);
  }    
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::insert(std::size_t entry_slot, std::size_t node_slot, bool first_insert) {
  return 0;
}

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
// RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::insert(const SlotIndT& entry_ind, const SlotIndT& node_ind, bool first_insert) {

//   // Current node and entry to add
//   TreeNode& node = m_nodes[node_ind];
//   TreeEntry& entry = m_entries[entry_ind];

//   // I4: adjust the bounding box of this node to include the newly added entry
//   node.extend(entry);

//   if(node.is_leaf) {

//     // CS2: if `node` is a leaf, insert the new `entry` into it
//     // (We take care of overfull leaves below)
//     node.add_child(entry_ind);
//   }
//   else {

//     // I1: call `choose_subtree` to find an appropriate node in which to place the entry, starting from the current `node_ind`, and insert it there
//     SlotIndT new_node_ind = insert(entry_ind,
// 				   choose_subtree(node_ind, entry_ind),
// 				   first_insert);

//     // No modifications need to be propagated up the tree
//     if(new_node_ind == NodePool::INVALID_SLOT) {     
//       return new_node_ind;
//     }

//     // The downstream insert resulted in a `new_node_ind` that now needs to be recorded
//     node.add_child(new_node_ind);
//   }

//   // If the procedure above resulted in an overfull tree node, handle it
//   if(node.num_child_nodes > MAX_NODESIZE) {
//     return overflow_treatment(node_ind, first_insert);
//   }
  
//   return NodePool::INVALID_SLOT;
// }

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
// RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::overflow_treatment(const SlotIndT& node_ind, bool first_insert) {
//   return NodePool::INVALID_SLOT;
// }

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
// RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::choose_subtree(const SlotIndT& start_node_ind, const SlotIndT& entry_ind) {
  
//   TreeNode& start_node = m_nodes[start_node_ind];
//   TreeEntry& entry = m_entries[entry_ind];

//   // This algorithm requires to start at an internal tree node (i.e. a node that does not directly contain entries)
//   assert(!start_node.is_leaf);

//   // Empty nodes cannot exist
//   assert(start_node.num_child_nodes > 0);

//   TreeNode& child = m_nodes[start_node.child_inds[0]];
//   if(child.is_leaf) {

//     // All child elements of `start_node` are leaves
//     // CS2: choose the entry whose bounding box needs the least overlap enlargement to accommodate the new data

//     // std::min_element()
//   }

//   // The child elements of `start_node` are internal nodes themselves
//   // CS2: choose the entry in `start_node` whose bounding box needs the least volume enlargement to include the new data
//   // (Recall that `start_node` already has its bounding box extended to include the new entry)

//   // std::min_element()
// }

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
// void RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::Rebuild() {

//   // implement STL algorithm to rebuild tree

//   // 1) get a list of all TreeEntries from the memory pool
//   // 2) reset both memory pools
//   // 3) sort the list into groups such that each group goes into one TreeEntry
//   // 4) rebuild the new tree structure one group at a time:
//   //    4.1) add the elements in one group to subsequent entries in the entry pool
//   //    4.2) insert the corresponding nodes into the node pool
  
//   // -> also need to reorder the memory pools  
// }

// // template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
// // const PayloadT& RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::Search(const IndexT& ind) {

// // }

