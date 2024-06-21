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
BoundingBox<CoordT, dims>::BoundingBox() : start_coords(std::numeric_limits<CoordT>::max()), end_coords(std::numeric_limits<CoordT>::min()),
					   shape((CoordT)0u) { };

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
void BoundingBox<CoordT, dims>::reset_bounding_box() {
  start_coords = std::numeric_limits<CoordT>::max();
  end_coords = std::numeric_limits<CoordT>::min();
  shape = (CoordT)0u;
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
Vector<CoordT, dims> BoundingBox<CoordT, dims>::compute_center() {
  return (start_coords + end_coords) / (CoordT)2;
}

template <class CoordT, std::size_t dims>
CoordT BoundingBox<CoordT, dims>::compute_center_distance(const BoundingBox<CoordT, dims>& bbox) {

  Vector<CoordT, dims> center_distance_vec = bbox.compute_center() - compute_center();

  // Actually calculates the squared distance, which is fine, as we only use it for relative comparisons
  return VectorUtils::inner(center_distance_vec, center_distance_vec);
}

template <class CoordT, std::size_t dims>
CoordT BoundingBox<CoordT, dims>::compute_center_distance(const Vector<CoordT, dims>& center) {

  Vector<CoordT, dims> center_distance_vec = center - compute_center();  
  return VectorUtils::inner(center_distance_vec, center_distance_vec);
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

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeEntry::TreeEntry() : BoundingBox<CoordT, dims>(), payload() { }

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::TreeNode::TreeNode() : BoundingBox<CoordT, dims>(), is_leaf(true), num_children(0) {
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
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::build_new_leaf_node(const std::vector<std::size_t>& entry_slots) {
  
  std::size_t node_slot = m_nodes.get_empty_slot_front();
  TreeNode& node = m_nodes[node_slot];

  // The new node is a leaf ...
  node.set_as_empty_leaf_node();

  // ... has a new bounding box ...
  node.reset_bounding_box();
  
  // ... and holds a reference to the new entry
  for(std::size_t entry_slot: entry_slots) {
    node.add_child(entry_slot);
  }
  
  return node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::build_new_internal_node(const std::vector<std::size_t>& child_node_slots) {
  
  std::size_t node_slot = m_nodes.get_empty_slot_front();
  TreeNode& node = m_nodes[node_slot];

  // The new node is an internal one ...
  node.set_as_empty_internal_node();

  // ... has a new bounding box ...
  node.reset_bounding_box();
  
  // ... and already comes with a few children
  for(std::size_t child_slot: child_node_slots) {
    node.add_child(child_slot);
  }
  
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
    m_root_slot = build_new_leaf_node({entry_slot});
  }
  else {

    // Non-empty tree: traverse the tree starting from the root node and insert the new entry at the correct location
    insert(entry_slot, m_root_slot);
  }    
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::insert(std::size_t entry_slot, std::size_t start_node_slot, bool first_insert) {  

  // Current node and entry to add
  TreeNode& start_node = m_nodes[start_node_slot];
  TreeEntry& entry = m_entries[entry_slot];
  
  // I4: adjust the bounding box of this node to include the newly added entry
  start_node.extend(entry);

  if(start_node.is_leaf) {

    // CS2: if `start_node` already is a leaf, simply insert the new `entry` into it
    // (We take care of overfull leaves below)
    start_node.add_child(entry_slot);
  }
  else {

    // I1: call `choose_subtree` to walk the tree and find the next `start_node` to check, and insert the `entry` there
    std::size_t new_node_slot = insert(entry_slot,
				       choose_subtree(entry_slot, start_node_slot),
				       first_insert);

    // Check if a new node was created during the insertion process
    if(new_node_slot == NodePool::INVALID_SLOT) {

      // No new slot created; nothing needs to be done by the caller
      return NodePool::INVALID_SLOT;
    }

    // The downstream insert resulted in a new node at `new_node_slot` that now needs to be added to our `start_node`
    start_node.add_child(new_node_slot);
  }

  // If the procedure above resulted in an overfull start_node, take care of it now
  if(start_node.num_children > MAX_NODESIZE) {
    return overflow_treatment(start_node_slot, first_insert);
  }
  
  return NodePool::INVALID_SLOT;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::choose_subtree(std::size_t entry_slot, std::size_t start_node_slot) {
  
  TreeNode& start_node = m_nodes[start_node_slot];
  TreeEntry& entry = m_entries[entry_slot];

  // Empty nodes cannot exist
  assert(start_node.num_children > 0);
  
  // This algorithm requires to start at an internal tree node (i.e. a node that does not directly contain entries)
  assert(!start_node.is_leaf);  

  TreeNode& child = m_nodes[start_node.child_slots[0]];
  if(child.is_leaf) {

    // All child nodes of `start_node` are leaf nodes
    // CS2: choose the entry in `start_node` whose bounding box needs the least overlap enlargement to accommodate the new data
    return find_min_overlap_enlargement_child(start_node_slot, entry);
  }

  // The child nodes of `start_node` are internal nodes themselves
  // CS2: choose the entry in `start_node` whose bounding box needs the least volume enlargement to include the new data
  return find_min_volume_enlargement_child(start_node_slot, entry);
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::find_min_overlap_enlargement_child(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add) {

  TreeNode& node = m_nodes[node_slot];

  // This only makes sense if the children of `node_slot` are themselves nodes that can accept `to_add`, i.e. `node_slot` is not a leaf
  assert(!node.is_leaf);

  // Returns `true` if adding `to_add` to the child at `child_slot_a` is strictly more advantageous than
  // adding it to `child_slot_b`.
  auto comp = [this, &node_slot, &to_add](const std::size_t& child_slot_a, const std::size_t& child_slot_b) -> bool {

    // Test for overlap enlargement first
    std::size_t overlap_enlargement_a = child_overlap_enlargement(node_slot, child_slot_a, to_add);
    std::size_t overlap_enlargement_b = child_overlap_enlargement(node_slot, child_slot_b, to_add);
    if(overlap_enlargement_a < overlap_enlargement_b) {
      return true;
    }
    else if(overlap_enlargement_a == overlap_enlargement_b) {

      // Break ties by looking at volume enlargement
      std::size_t volume_enlargement_a = node_volume_enlargement(child_slot_a, to_add);
      std::size_t volume_enlargement_b = node_volume_enlargement(child_slot_b, to_add);
      
      return volume_enlargement_a < volume_enlargement_b;
    }
    
    return false;    
  };

  return *std::min_element(node.child_slots.begin(),
			   node.child_slots.begin() + node.num_children,
			   comp); 
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::find_min_volume_enlargement_child(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add) {

  TreeNode& node = m_nodes[node_slot];
  
  // This only makes sense if the children of `node_slot` are themselves nodes that can accept `to_add`, i.e. `node_slot` is not a leaf
  assert(!node.is_leaf);

  // Returns `true` if adding `to_add` to the child at `child_slot_a` is strictly more advantagous than
  // adding it to `child_slot_b`
  auto comp = [this, &to_add](const std::size_t& child_slot_a, const std::size_t& child_slot_b) -> bool {

    // Test for volume enlargement first
    std::size_t volume_enlargement_a = node_volume_enlargement(child_slot_a, to_add);
    std::size_t volume_enlargement_b = node_volume_enlargement(child_slot_b, to_add);
    if(volume_enlargement_a < volume_enlargement_b) {
      return true;
    }
    else if(volume_enlargement_a == volume_enlargement_b) {

      // Break ties by looking at the absolute volume
      return node_volume(child_slot_a) < node_volume(child_slot_b);
    }
    
    return false;      
  };

  return *std::min_element(node.child_slots.begin(),
			   node.child_slots.begin() + node.num_children,
			   comp);
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::child_overlap_enlargement(std::size_t node_slot, std::size_t child_slot,
												 const BoundingBox<CoordT, dims>& to_add) {
  // Get the bounding box of the child node ...
  BoundingBox<CoordT, dims>& child_bbox = get_bbox(node_slot, child_slot);

  // ... and also the bounding box of the child, extended by `to_add`
  BoundingBox<CoordT, dims> extended_child_bbox(child_bbox);
  extended_child_bbox.extend(to_add);

  TreeNode& node = m_nodes[node_slot];
  
  std::size_t child_overlap_enlargement = 0u;  
  for(std::size_t i = 0; i < node.num_children; i++) {
    
    std::size_t other_child_slot = node.child_slots[i];
    if(other_child_slot == child_slot) {
      continue;  // Don't compute overlap against itself
    }
    
    BoundingBox<CoordT, dims>& other_child_bbox = get_bbox(node_slot, other_child_slot);
    child_overlap_enlargement += (extended_child_bbox.compute_overlap_volume(other_child_bbox) - child_bbox.compute_overlap_volume(other_child_bbox));
  }

  return child_overlap_enlargement;  
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::node_volume_enlargement(std::size_t node_slot, const BoundingBox<CoordT, dims>& to_add) {

  BoundingBox<CoordT, dims>& node_bbox = m_nodes[node_slot];
  BoundingBox<CoordT, dims> extended_node_bbox(node_bbox);
  extended_node_bbox.extend(to_add);

  return extended_node_bbox.volume() - node_bbox.volume();
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::node_volume(std::size_t node_slot) {
  return m_nodes[node_slot].volume();
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
BoundingBox<CoordT, dims>& RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::get_bbox(std::size_t node_slot, std::size_t child_slot) {
  
  TreeNode& node = m_nodes[node_slot];
  if(node.is_leaf) {
    
    // This node's children are actually entries
    return m_entries[child_slot];
  }
  else {

    // This node's children are other nodes
    return m_nodes[child_slot];
  }
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::recalculate_bbox(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];
  node.reset_bounding_box(); // reset everything ...

  // ... and grow the bounding box again for all children
  for(std::size_t i = 0; i < node.num_children; i++) {    
    std::size_t cur_child_slot = node.child_slots[i];    
    node.extend(get_bbox(node_slot, cur_child_slot));
  }
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::overflow_treatment(std::size_t node_slot, bool first_insert) {

  // This assumes we are given a node that is overfull and needs to be cleaned up
  assert(m_nodes[node_slot].num_children == MAX_NODESIZE + 1);
  
  // OT1: If we are sitting at a leaf node that is not the root _and_ this is the first call during an insertion, use reinsertion to clean up the overfull node ...
  // (Note: this is a deviation from the original R*-Tree, which permits forced reinsertion at all levels of the tree, not just for root nodes)
  if((node_slot != m_root_slot) && m_nodes[node_slot].is_leaf && first_insert) {
    
    reinsert(node_slot);
    return NodePool::INVALID_SLOT;
  }

  // OT1: ... otherwise split the node to distribute its too many children
  std::size_t new_node_slot = split(node_slot);

  // If we actually split the root node, create a new root node whose children are `node_slot` and `new_node_slot`
  if(node_slot == m_root_slot) {

    std::size_t new_root_slot = build_new_internal_node({node_slot, new_node_slot});

    // I4: recalculate the bounding box of the new root node
    recalculate_bbox(new_root_slot);

    // Update the root node
    m_root_slot = new_root_slot;
    
    // Already updated the root node, nothing to propagate upwards
    return NodePool::INVALID_SLOT;
  }

  return new_node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::reinsert(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];

  // We have a leaf node with too many entries and are trying to re-insert some of them into the tree
  assert(node.is_leaf);
  
  // We have an overfull node that we need to clean up
  assert(node.num_children == MAX_NODESIZE + 1);

  constexpr std::size_t num_children = MAX_NODESIZE + 1;
  constexpr std::size_t p = std::max<std::size_t>(REINSERT_P_FRAC * num_children, 1u);
  assert(p < num_children);

  // RI1 + RI2: Sort the children of the node at `node_slot` in increasing order w.r.t. the distance between
  // the center of the child's bounding rectangle and the center of the bounding rectangle of the node
  Vector<CoordT, dims> node_center = node.compute_center();
  auto comp = [this, &node_slot, &node_center](const std::size_t& child_slot_a, const std::size_t& child_slot_b) -> bool {

    CoordT center_distance_a = get_bbox(node_slot, child_slot_a).compute_center_distance(node_center);
    CoordT center_distance_b = get_bbox(node_slot, child_slot_b).compute_center_distance(node_center);
    
    return center_distance_a < center_distance_b;
  };
  std::partial_sort(node.child_slots.begin(), node.child_slots.end() - p, node.child_slots.end(), comp);

  // After the partial sort, the first `num_children - p` elements contain the `num_children - p` entries whose
  // bounding box centers have the smallest distance to the center of the bounding box of the node  
  
  // Remove the `p` children with the largest distance from this node ...
  std::array<std::size_t, p> removed_entry_node_slots;
  std::copy(node.child_slots.end() - p, node.child_slots.end(), removed_entry_node_slots.begin());
  std::fill(node.child_slots.end() - p, node.child_slots.end(), NodePool::INVALID_SLOT);
  node.num_children = num_children - p;
  recalculate_bbox(node_slot);

  // ... and reinsert them into the tree starting at the root node
  for(std::size_t entry_slot_to_insert: removed_entry_node_slots) {
    insert(entry_slot_to_insert, m_root_slot, false);
  }
}

template <typename CoordT, std::size_t dims, class PayloadT, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
std::size_t RTree<CoordT, dims, PayloadT, MAX_NODESIZE, MIN_NODESIZE>::split(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];

  // We have an overfull node that we need to clean up
  assert(node.num_children == MAX_NODESIZE + 1);
  
  // Determine how to best split the entries into two groups, each of which contains less than `MAX_NODESIZE` elements

  // First, create a new node that will take over some of the elements
  // std::size_t new_node_slot = 

  
  
  return NodePool::INVALID_SLOT;
}

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

