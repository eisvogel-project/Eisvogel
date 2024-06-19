#include <algorithm>

template <class T, std::unsigned_integral SlotIndT>
MemoryPool<T, SlotIndT>::MemoryPool(std::size_t init_size) : m_init_size(init_size),
							 m_free_start(Slot::INVALID), m_free_end(Slot::INVALID),
							 m_alloc_start(Slot::INVALID), m_alloc_end(Slot::INVALID),
							 m_data(init_size), m_slots(init_size) {
  assert(m_init_size > 0);
  reset_slot_lists();
}

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::reset() {

  m_data.clear();
  m_data.resize(m_init_size);

  m_slots.clear();
  m_slots.resize(m_init_size);
  
  reset_slot_lists();
}

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::reset_slot_lists() {

  assert(m_data.size() == m_slots.size());

  // Reset the two doubly-linked lists
  m_alloc_start = Slot::INVALID;
  m_alloc_end = Slot::INVALID;

  m_free_start = Slot::INVALID;
  m_free_end = Slot::INVALID;
  
  // Mark everything as free
  for(SlotIndT i = 0; i < m_slots.size(); i++) {
    add_to_slot_list_back(i, m_free_start, m_free_end);
  }
}

template <class T, std::unsigned_integral SlotIndT>
T& MemoryPool<T, SlotIndT>::operator[](SlotIndT ind) {
  return m_data[ind];
}

template <class T, std::unsigned_integral SlotIndT>
const T& MemoryPool<T, SlotIndT>::operator[](SlotIndT ind) const {
  return m_data[ind];
}

template <class T, std::unsigned_integral SlotIndT>
SlotIndT MemoryPool<T, SlotIndT>::get_empty_slot_front() {

  // Grow the pool if needed
  if(m_free_start == Slot::INVALID) {
    grow();
  }

  // Get the index of the next free slot
  SlotIndT alloc_slot_ind = m_free_start;  
  
  // Remove the newly allocated slot from the beginning of the list of free slots ...
  remove_from_slot_list(alloc_slot_ind, m_free_start, m_free_end);

  // ... and add it to the back of the list of allocated slots
  add_to_slot_list_back(alloc_slot_ind, m_alloc_start, m_alloc_end);

  // Return the allocated slot
  return alloc_slot_ind;
}

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::free_slot(SlotIndT ind) {

  // Remove the slot with index `ind` from the list of allocated slots ...
  remove_from_slot_list(ind, m_alloc_start, m_alloc_end);

  // ... and add it to the back of the list of free slots
  add_to_slot_list_back(ind, m_free_start, m_free_end);
}

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::grow() {

  // Double the size of the currently-allocated memory
  std::size_t current_size = m_data.size();
  std::size_t new_size = 2 * current_size;
  m_data.resize(new_size);

  // Extend the slot list
  std::vector<Slot> new_slots_free(new_size - current_size);
  m_slots.insert(m_slots.end(), new_slots_free.begin(), new_slots_free.end());

  // Mark the new slots as free
  for(SlotIndT i = current_size; i < new_size; i++) {
    add_to_slot_list_back(i, m_free_start, m_free_end);
  }  
}

template <class T, std::unsigned_integral SlotIndT>
bool MemoryPool<T, SlotIndT>::is_valid_slot(SlotIndT& ind) {
  return (ind != Slot::INVALID) && (ind < m_slots.size());
}

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::remove_from_slot_list(SlotIndT& to_remove, SlotIndT& list_start, SlotIndT& list_end) {

  assert(is_valid_slot(to_remove));
  
  // Get the indices of the previous and next element from the list
  SlotIndT& next_elem = m_slots[to_remove].next_data_ind;
  SlotIndT& prev_elem = m_slots[to_remove].prev_data_ind;

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

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::add_to_slot_list_back(SlotIndT& to_add, SlotIndT& list_start, SlotIndT& list_end) {

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

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::add_to_slot_list_front(SlotIndT& to_add, SlotIndT& list_start, SlotIndT& list_end) {

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

template <class T, std::unsigned_integral SlotIndT>
template <class CallableT>
void MemoryPool<T, SlotIndT>::loop_over_elements_front_to_back(SlotIndT& list_start, CallableT&& worker) {

  SlotIndT cur_elem = list_start;
  while(cur_elem != Slot::INVALID) {
    
    // call worker
    worker(cur_elem);
    // worker(m_data[cur_elem]);

    // Walk to the following element in the list
    cur_elem = m_slots[cur_elem].next_data_ind;
  }  
}

template <class T, std::unsigned_integral SlotIndT>
template <class CallableT>
void MemoryPool<T, SlotIndT>::loop_over_elements_back_to_front(SlotIndT& list_end, CallableT&& worker) {

  SlotIndT cur_elem = list_end;
  while(cur_elem != Slot::INVALID) {

    // call worker
    worker(cur_elem);
    // worker(m_data[cur_elem]);

    // Walk to the previoud element in the list
    cur_elem = m_slots[cur_elem].prev_data_ind;
  }
}

// ------

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::print_free() {

  auto printer = [](const SlotIndT& cur_slot) {
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

template <class T, std::unsigned_integral SlotIndT>
void MemoryPool<T, SlotIndT>::print_allocated() {

  auto printer = [](const SlotIndT& cur_slot) {
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

template <class T, std::unsigned_integral SlotIndT>
bool MemoryPool<T, SlotIndT>::is_free(SlotIndT ind) {
  return false; // TODO
}

template <class T, std::unsigned_integral SlotIndT>
bool MemoryPool<T, SlotIndT>::is_allocated(SlotIndT ind) {
  return !is_free(ind);
}

// =======================

template <class IndexT, std::size_t dims>
bool BoundingBox<IndexT, dims>::contains(const IndexT& ind) {
  return false;
}

template <class IndexT, std::size_t dims>
bool BoundingBox<IndexT, dims>::overlaps(const BoundingBox& bbox) {
  return false;
}

template <class IndexT, std::size_t dims>
void BoundingBox<IndexT, dims>::stretch(const BoundingBox& bbox) {

}

template <class IndexT, std::size_t dims>
std::size_t BoundingBox<IndexT, dims>::compute_overlap(const BoundingBox& other) {
  return 0;
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

// Default constructor: mark everything as invalid
template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::Node() : is_leaf(true), num_child_nodes(0) {
  child_inds.fill(0);
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::mark_as_empty_leaf() {
  is_leaf = true;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::mark_as_empty_internal() {
  is_leaf = false;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::add_child(SlotIndT child_ind) {
  assert(num_child_nodes <= MAX_NODESIZE);
  child_inds[num_child_nodes] = child_ind;
  num_child_nodes++;
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
SlotIndT Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::get_min_area_enlargement_child(SlotIndT& new_entry_ind) {

  
  for(std::size_t i = 0; i < num_child_nodes; i++) {
    
  }
  
  return 0;
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
std::size_t Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::get_child_overlap(SlotIndT& child_ind) {

  return 0;
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
SlotIndT Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::get_min_overlap_enlargement_child(SlotIndT& new_entry) {

  SlotIndT min_child = std::nullopt;
  
  return min_child;
}

template <class IndexT, std::size_t dims, typename SlotIndT, std::size_t MAX_NODESIZE>
template <class ComparatorT>
SlotIndT Node<IndexT, dims, SlotIndT, MAX_NODESIZE>::min_child(ComparatorT&& comp) {

  
  
  return 0;
};										 

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::RTree(std::size_t init_slot_size) : m_root_node_ind(NodePool::INVALID_SLOT),
											       m_nodes(init_slot_size), m_entries(init_slot_size) { }

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::InsertElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind) {

  // Take ownership of the `elem` and store it as an entry in the memory pool
  SlotIndT entry_ind = m_entries.get_empty_slot_front();
  TreeEntry& entry = m_entries[entry_ind];
  entry.payload = elem;
  entry.start_ind = start_ind;
  entry.end_ind = end_ind;

  // Insert the entry into the node structure of the tree
  if(m_root_node_ind == NodePool::INVALID_SLOT) {

    // Empty tree, create a new root node
    SlotIndT node_ind = m_nodes.get_empty_slot_front();
    m_root_node_ind = node_ind;
    TreeNode& node = m_nodes[node_ind];

    // The new root node is a leaf ...
    node.mark_as_empty_leaf();

    // ... and holds a reference to the new entry
    node.add_child(entry_ind);
  }
  else {

    // Non-empty tree: traverse the tree starting from the root node and insert the new entry at the correct location
    insert(entry_ind, m_root_node_ind);
  }    
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::insert(const SlotIndT& entry_ind, const SlotIndT& node_ind, bool first_insert) {

  // Current node and entry to add
  TreeNode& node = m_nodes[node_ind];
  TreeEntry& entry = m_entries[entry_ind];

  // I4: adjust the bounding box of this node to include the newly added entry
  node.stretch(entry);

  if(node.is_leaf) {

    // CS2: if `node` is a leaf, insert the new `entry` into it
    // (We take care of overfull leaves below)
    node.add_child(entry_ind);
  }
  else {

    // I1: call `choose_subtree` to find an appropriate node in which to place the entry, starting from the current `node_ind`, and insert it there
    SlotIndT new_node_ind = insert(entry_ind,
				   choose_subtree(node_ind, entry_ind),
				   first_insert);

    // No modifications need to be propagated up the tree
    if(new_node_ind == NodePool::INVALID_SLOT) {     
      return new_node_ind;
    }

    // The downstream insert resulted in a `new_node_ind` that now needs to be recorded
    node.add_child(new_node_ind);
  }

  // If the procedure above resulted in an overfull tree node, handle it
  if(node.num_child_nodes > MAX_NODESIZE) {
    return overflow_treatment(node_ind, first_insert);
  }
  
  return NodePool::INVALID_SLOT;
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::overflow_treatment(const SlotIndT& node_ind, bool first_insert) {
  return NodePool::INVALID_SLOT;
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::SlotIndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::choose_subtree(const SlotIndT& start_node_ind, const SlotIndT& entry_ind) {
  
  TreeNode& start_node = m_nodes[start_node_ind];
  TreeEntry& entry = m_entries[entry_ind];

  // This algorithm requires to start at an internal tree node (i.e. a node that does not directly contain entries)
  assert(!start_node.is_leaf);

  // Empty nodes cannot exist
  assert(start_node.num_child_nodes > 0);

  TreeNode& child = m_nodes[start_node.child_inds[0]];
  if(child.is_leaf) {

    // All child elements of `start_node` are leaves
    // CS2: choose the entry whose bounding box needs the least overlap enlargement to accommodate the new data

    // std::min_element()
  }

  // The child elements of `start_node` are internal nodes themselves
  // CS2: choose the entry in `start_node` whose bounding box needs the least volume enlargement to include the new data
  // (Recall that `start_node` already has its bounding box extended to include the new entry)

  // std::min_element()
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::Rebuild() {

  // implement STL algorithm to rebuild tree

  // 1) get a list of all TreeEntries from the memory pool
  // 2) reset both memory pools
  // 3) sort the list into groups such that each group goes into one TreeEntry
  // 4) rebuild the new tree structure one group at a time:
  //    4.1) add the elements in one group to subsequent entries in the entry pool
  //    4.2) insert the corresponding nodes into the node pool
  
  // -> also need to reorder the memory pools  
}

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
// const PayloadT& RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::Search(const IndexT& ind) {

// }

