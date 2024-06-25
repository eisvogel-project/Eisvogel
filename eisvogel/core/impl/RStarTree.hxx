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
bool BoundingBox<CoordT, dims>::overlaps(const BoundingBox<CoordT, dims>& bbox) const {

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
CoordT BoundingBox<CoordT, dims>::margin() {
  return std::accumulate(shape.cbegin(), shape.cend(), (CoordT)0u, std::plus<CoordT>());
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
  stream << "shape: " << bbox.shape     << "\n";
  stream << "-----------------------------------------\n";

  return stream;
}

// =======================

template <typename CoordT, std::size_t dims, class PayloadT>
RStarTree<CoordT, dims, PayloadT>::TreeEntry::TreeEntry() : ElemBoundingBox(), payload() { }

template <typename CoordT, std::size_t dims, class PayloadT>
RStarTree<CoordT, dims, PayloadT>::TreeNode::TreeNode() : ElemBoundingBox(), level(0), num_children(0) {
  child_slots.fill(0);
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::TreeNode::add_child(std::size_t child_slot) {
  assert(num_children <= MAX_NODESIZE);
  child_slots[num_children] = child_slot;
  num_children++;
}

template <typename CoordT, std::size_t dims, class PayloadT>
RStarTree<CoordT, dims, PayloadT>::RStarTree(std::size_t init_slot_size) : m_root_slot(MemoryPool<TreeNode>::INVALID_SLOT),
								   m_nodes(init_slot_size), m_entries(init_slot_size) { }

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::build_new_entry(const PayloadT& elem,
							   const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) {
  std::size_t entry_slot = m_entries.get_empty_slot_front();
  TreeEntry& entry = m_entries[entry_slot];
  
  entry.payload = elem;
  entry.start_coords = start_coords;
  entry.end_coords = end_coords;
  entry.shape = end_coords - start_coords;
  
  return entry_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::build_new_node(std::size_t level, const std::vector<std::size_t>& entry_slots) {

  std::size_t node_slot = m_nodes.get_empty_slot_front();
  TreeNode& node = m_nodes[node_slot];

  // Set the node attributes ... 
  node.reset_bounding_box();
  node.level = level;

  // ... and add the children
  node.num_children = 0;
  for(std::size_t entry_slot: entry_slots) {
    node.add_child(entry_slot);
  }

  return node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::InsertElement(const PayloadT& elem,
						  const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) {

  // Take ownership of the `elem` and store it as an entry in the memory pool
  std::size_t entry_slot = build_new_entry(elem, start_coords, end_coords);

  // Insert the entry into the node structure of the tree
  if(m_root_slot == NodePool::INVALID_SLOT) {
    
    // Empty tree: build new root node
    m_root_slot = build_new_node(0, {entry_slot});
    recalculate_bbox(m_root_slot);
  }
  else {

    // Non-empty tree: traverse the tree starting from the root node and insert the new entry at the correct location at level 0 (i.e. as a leaf node)
    insert_slot(entry_slot, 0, m_root_slot, true);
  }
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::insert_slot(std::size_t slot_to_add, std::size_t level_to_insert,
						       std::size_t start_node_slot, bool first_insert) {  

  // Start node: must be above the target level in the tree hierarchy!
  TreeNode& start_node = m_nodes[start_node_slot];  
  assert(start_node.level >= level_to_insert);
  
  ElemBoundingBox& bbox_to_add = get_bbox(slot_to_add, level_to_insert);
  
  // I4: adjust the bounding box of this node to include the newly added entry
  start_node.extend(bbox_to_add);

  if(start_node.level == level_to_insert) {

    // CS2: This node is already at the correct level, add the new entry to its children
    // (We take care of overfull leaves below)
    start_node.add_child(slot_to_add);
  }
  else {

    // I1: call `choose_subtree` to walk the tree and find the next `start_node` to check, and insert the `entry` there
    std::size_t new_node_slot = insert_slot(slot_to_add,
					    level_to_insert,
					    choose_subtree(slot_to_add, level_to_insert, start_node_slot),
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

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::choose_subtree(std::size_t slot_to_add, std::size_t level_to_insert, std::size_t start_node_slot) {
  
  TreeNode& start_node = m_nodes[start_node_slot];
  ElemBoundingBox& bbox_to_add = get_bbox(slot_to_add, level_to_insert);

  // Empty nodes cannot exist
  assert(start_node.num_children > 0);
  
  // This algorithm requires to start at an internal tree node (i.e. a node that does not directly contain entries)
  assert(start_node.level > 0);

  TreeNode& child = m_nodes[start_node.child_slots[0]];
  if(child.level == 0) {

    // All child nodes of `start_node` are leaf nodes
    // CS2: choose the entry in `start_node` whose bounding box needs the least overlap enlargement to accommodate the new data
    return find_min_overlap_enlargement_child(start_node_slot, bbox_to_add);
  }

  // The child nodes of `start_node` are internal nodes themselves
  // CS2: choose the entry in `start_node` whose bounding box needs the least volume enlargement to include the new data
  return find_min_volume_enlargement_child(start_node_slot, bbox_to_add);
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::find_min_overlap_enlargement_child(std::size_t node_slot, const ElemBoundingBox& to_add) {

  TreeNode& node = m_nodes[node_slot];

  // This only makes sense if the children of `node_slot` are themselves nodes that can accept `to_add`, i.e. `node_slot` is not a leaf
  assert(node.level > 0);

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

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::find_min_volume_enlargement_child(std::size_t node_slot, const ElemBoundingBox& to_add) {

  TreeNode& node = m_nodes[node_slot];
  
  // This only makes sense if the children of `node_slot` are themselves nodes that can accept `to_add`, i.e. `node_slot` is not a leaf
  assert(node.level > 0);

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

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::child_overlap_enlargement(std::size_t node_slot, std::size_t child_slot,
												 const ElemBoundingBox& to_add) {
  TreeNode& node = m_nodes[node_slot];
    
  // Get the bounding box of the child node ...
  ElemBoundingBox& child_bbox = get_bbox(child_slot, node.level);

  // ... and also the bounding box of the child, extended by `to_add`
  ElemBoundingBox extended_child_bbox(child_bbox);
  extended_child_bbox.extend(to_add);
  
  std::size_t child_overlap_enlargement = 0u;  
  for(std::size_t i = 0; i < node.num_children; i++) {
    
    std::size_t other_child_slot = node.child_slots[i];
    if(other_child_slot == child_slot) {
      continue;  // Don't compute overlap against itself
    }
    
    ElemBoundingBox& other_child_bbox = get_bbox(other_child_slot, node.level);
    child_overlap_enlargement += (extended_child_bbox.compute_overlap_volume(other_child_bbox) - child_bbox.compute_overlap_volume(other_child_bbox));
  }

  return child_overlap_enlargement;  
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::node_volume_enlargement(std::size_t node_slot, const ElemBoundingBox& to_add) {

  ElemBoundingBox& node_bbox = m_nodes[node_slot];
  ElemBoundingBox extended_node_bbox(node_bbox);
  extended_node_bbox.extend(to_add);

  return extended_node_bbox.volume() - node_bbox.volume();
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::node_volume(std::size_t node_slot) {
  return m_nodes[node_slot].volume();
}

template <typename CoordT, std::size_t dims, class PayloadT>
RStarTree<CoordT, dims, PayloadT>::ElemBoundingBox&
RStarTree<CoordT, dims, PayloadT>::get_bbox(std::size_t slot, std::size_t level) {

  if(level == 0) {
    return m_entries[slot];
  }
  else {
    return m_nodes[slot];
  }    
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::recalculate_bbox(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];
  node.reset_bounding_box(); // reset everything ...

  // ... and grow the bounding box again for all children
  for(std::size_t i = 0; i < node.num_children; i++) {    
    std::size_t cur_child_slot = node.child_slots[i];    
    node.extend(get_bbox(cur_child_slot, node.level));
  }
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::overflow_treatment(std::size_t node_slot, bool first_insert) {

  std::cout << "-------" << std::endl;
  std::cout << "in overflow_treatment" << std::endl;
  
  // This assumes we are given a node that is overfull and needs to be cleaned up
  assert(m_nodes[node_slot].num_children == MAX_NODESIZE + 1);
  
  // OT1: If we are sitting at a node that is not the root _and_ this is the first call during an insertion, use reinsertion at this level to clean up the overfull node ...
  if((node_slot != m_root_slot) && first_insert) {

    std::cout << "treat overflow by reinsertion" << std::endl;
    
    reinsert(node_slot);
    return NodePool::INVALID_SLOT;
  }

  std::cout << "treat overflow by split" << std::endl;
  
  // OT1: ... otherwise split the node to distribute its too many children
  std::size_t new_node_slot = split(node_slot);

  std::cout << "split gave new node at slot = " << new_node_slot << std::endl;

  std::cout << "situation after split:" << std::endl;
  std::cout << "node_slot:" << std::endl;
  std::cout << m_nodes[node_slot] << std::endl;

  std::cout << "new_node_slot:" << std::endl;
  std::cout << m_nodes[new_node_slot] << std::endl;
  
  // If we actually split the root node, create a new root node whose children are `node_slot` and `new_node_slot`
  if(node_slot == m_root_slot) {

    std::size_t current_root_level = m_nodes[m_root_slot].level;
    std::size_t new_root_slot = build_new_node(current_root_level + 1, {node_slot, new_node_slot});

    // I4: recalculate the bounding box of the new root node
    recalculate_bbox(new_root_slot);

    // Update the root node
    m_root_slot = new_root_slot;

    std::cout << "new_root_slot:" << std::endl;
    std::cout << m_nodes[new_root_slot] << std::endl;
    
    // Already updated the root node, nothing to propagate upwards
    return NodePool::INVALID_SLOT;
  }

  std::cout << "-------" << std::endl;
  
  return new_node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::reinsert(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];
  
  // We have an overfull node that we need to clean up
  assert(node.num_children == MAX_NODESIZE + 1);

  constexpr std::size_t num_children = MAX_NODESIZE + 1;
  constexpr std::size_t p = std::max<std::size_t>(REINSERT_P_FRAC * num_children, 1u);
  assert(p < num_children);

  // RI1 + RI2: Sort the children of the node at `node_slot` in increasing order w.r.t. the distance between
  // the center of the child's bounding rectangle and the center of the bounding rectangle of the node
  Vector<CoordT, dims> node_center = node.compute_center();
  auto comp = [this, &node, &node_center](const std::size_t& child_slot_a, const std::size_t& child_slot_b) -> bool {

    CoordT center_distance_a = get_bbox(child_slot_a, node.level).compute_center_distance(node_center);
    CoordT center_distance_b = get_bbox(child_slot_b, node.level).compute_center_distance(node_center);
    
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

  // ... and reinsert them into the tree at the same level from which they were removed, starting at the root node
  for(std::size_t entry_slot_to_insert: removed_entry_node_slots) {
    insert_slot(entry_slot_to_insert, node.level, m_root_slot, false);
  }
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::size_t RStarTree<CoordT, dims, PayloadT>::split(std::size_t node_slot) {

  TreeNode& node = m_nodes[node_slot];

  // We have an overfull node that we need to clean up
  assert(node.num_children == MAX_NODESIZE + 1);
  
  // Number of possible ways in which to split the children into two groups, each of which contains less than `MAX_NODESIZE` elements,
  // but at least `MIN_NODESIZE` elements
  constexpr std::size_t number_possible_splits = MAX_NODESIZE + 2 - 2 * MIN_NODESIZE;

  // Bounding boxes to use when testing various splits 
  ElemBoundingBox bbox_1, bbox_2;

  enum Edge {
    Lower, Upper
  };
  
  // Parameterization of the best way to split the entries in this node
  std::size_t best_split_axis = 0;
  Edge best_edge = Edge::Lower;
  std::size_t best_split = 0;
  
  CoordT best_split_axis_margin = std::numeric_limits<CoordT>::max();
  
  // The following is a combined determination of the axis along which the node split is best performed, and the
  // partitioning of the node entries that is best
  for(std::size_t cur_axis = 0; cur_axis < dims; cur_axis++) {
    
    // The margin for this particular axis (calculated below) is used to identify the axis along which to split
    CoordT cur_axis_margin = (CoordT)0;

    // Also keep track of the overlap- and volume-values for a particular split
    CoordT cur_axis_min_overlap = std::numeric_limits<CoordT>::max();
    CoordT cur_axis_min_volume = std::numeric_limits<CoordT>::max();
    Edge cur_axis_best_edge = Edge::Lower;
    std::size_t cur_axis_best_split = 0;
    
    // Identifying the best split involves iterating over different ways of sorting the entries ... 
    for(auto cur_edge: {Edge::Lower, Edge::Upper}) {
        
        if(cur_edge == Edge::Lower) {

	  // Sort nodes by the lower values of their bounding boxes along `cur_axis`
	  auto sort_lower = [this, &node, &cur_axis](const std::size_t& child_slot_a, const std::size_t& child_slot_b) {
	    return get_bbox(child_slot_a, node.level).start_coords[cur_axis] < get_bbox(child_slot_b, node.level).start_coords[cur_axis];
	  };
	  std::sort(node.child_slots.begin(), node.child_slots.end(), sort_lower);
        }
	else {
	  
	  // Sort nodes by the upper values of their bounding boxes along `cur_axis`
	  auto sort_upper = [this, &node, &cur_axis](const std::size_t& child_slot_a, const std::size_t& child_slot_b) {
	    return get_bbox(child_slot_a, node.level).end_coords[cur_axis] < get_bbox(child_slot_b, node.level).end_coords[cur_axis];
	  };
	  std::sort(node.child_slots.begin(), node.child_slots.end(), sort_upper);
        }

	// ... and different ways of actually splitting the entries
	for(std::size_t cur_split = 0; cur_split < number_possible_splits; cur_split++) {

	  // The split we're currently investigating is [0, MIN_NODESIZE + cur_split] and [MIN_NODESIZE + cur_split + 1, MAX_NODESIZE + 1]
	  // Calculate the bounding boxes of the two partitions
	  bbox_1.reset_bounding_box();
	  auto extend_bbox_1 = [this, &node, &bbox_1](const std::size_t& child_slot) {bbox_1.extend(get_bbox(child_slot, node.level));};
	  std::for_each(node.child_slots.begin(), node.child_slots.begin() + MIN_NODESIZE + cur_split, extend_bbox_1);
	  
	  auto extend_bbox_2 = [this, &node, &bbox_2](const std::size_t& child_slot) {bbox_2.extend(get_bbox(child_slot, node.level));};
	  bbox_2.reset_bounding_box();
	  std::for_each(node.child_slots.begin() + MIN_NODESIZE + cur_split, node.child_slots.end(), extend_bbox_2);

	  // Compute statistics for this split
	  cur_axis_margin += (bbox_1.margin() + bbox_2.margin());
	  CoordT cur_split_volume = bbox_1.volume() + bbox_2.volume();
	  CoordT cur_split_overlap = bbox_1.compute_overlap_volume(bbox_2);

	  // Keep updating the best split for this axis
	  if((cur_split_overlap < cur_axis_min_overlap) ||
	     ((cur_split_overlap == cur_axis_min_overlap) && (cur_split_volume < cur_axis_min_volume))) {
	    
	    cur_axis_min_overlap = cur_split_overlap;
	    cur_axis_min_volume = cur_split_volume;

	    cur_axis_best_split = cur_split;
	    cur_axis_best_edge = cur_edge;
	  }
	}	
    }

    // Keep track of the axis which minimum margin; this is the one along which the split will be performed below
    if(cur_axis_margin < best_split_axis_margin) {
      best_split_axis_margin = cur_axis_margin;
      best_split_axis = cur_axis;
      best_edge = cur_axis_best_edge;
      best_split = cur_axis_best_split;
    }
  }

  // Now, all splits have been examined, and the best split axis, edge for the sort-before-split, and partitioning have been determined

  // Reconstruct the order that was found to be the optimal one
  if(best_edge == Edge::Lower) {

    // Re-sort according to the lower edge along the split axis
    auto sort_lower = [this, &node, &best_split_axis](const std::size_t& child_slot_a, const std::size_t& child_slot_b) {
      return get_bbox(child_slot_a, node.level).start_coords[best_split_axis] < get_bbox(child_slot_b, node.level).start_coords[best_split_axis];
    };
    std::sort(node.child_slots.begin(), node.child_slots.end(), sort_lower);    
  }
  else if(best_split_axis != dims - 1) {

    // Re-sort according to the upper edge along the split axis
    // (But don't need to re-sort if the state after the above loop is already the best, in which case everything is already correct)
    auto sort_upper = [this, &node, &best_split_axis](const std::size_t& child_slot_a, const std::size_t& child_slot_b) {
      return get_bbox(child_slot_a, node.level).end_coords[best_split_axis] < get_bbox(child_slot_b, node.level).end_coords[best_split_axis];
    };
    std::sort(node.child_slots.begin(), node.child_slots.end(), sort_upper);
  }

  // The state of the node and its children is now in the correct state for the best split

  // Extract the children that will go into the new node ...
  std::vector<std::size_t> new_node_children(node.child_slots.begin() + MIN_NODESIZE + best_split, node.child_slots.end());

  // ... and remove them from the old node
  std::fill(node.child_slots.begin() + MIN_NODESIZE + best_split, node.child_slots.end(), NodePool::INVALID_SLOT);
  node.num_children = MIN_NODESIZE + best_split;
  
  // Create a new node at the same tree level which will take over some of the excess elements
  std::size_t new_node_slot = build_new_node(node.level, new_node_children);

  // Recalculate the bounding boxes of both the old and the new node
  recalculate_bbox(new_node_slot);
  recalculate_bbox(node_slot);
  
  return new_node_slot;
}

template <typename CoordT, std::size_t dims, class PayloadT>
const PayloadT& RStarTree<CoordT, dims, PayloadT>::Search(const Vector<CoordT, dims>& coords) {
  
  std::size_t cur_node_slot = m_root_slot;
  while(true) {

    TreeNode& cur_node = m_nodes[cur_node_slot];
    
    // Iterate over all children in this node
    for(std::size_t i = 0; i < cur_node.num_children; i++) {
      std::size_t cur_child_slot = cur_node.child_slots[i];
            
      // This is the correct child to follow
      if(get_bbox(cur_child_slot, cur_node.level).contains(coords)) {
	
	if(cur_node.level == 0) {
	  
	  // Already at a leaf node, found the entry!
	  return m_entries[cur_child_slot].payload;
	}
	else {

	  // Not yet at the leaf level; continue with the child
	  cur_node_slot = cur_child_slot;
	  break;
	}
      }
    }

    if(cur_node.level == 0) {
      throw std::runtime_error("Element not found!");
    }
  } 
}

template <typename CoordT, std::size_t dims, class PayloadT>
std::vector<std::reference_wrapper<const PayloadT>> RStarTree<CoordT, dims, PayloadT>::Search(const Vector<CoordT, dims>& start_coords, const Vector<CoordT, dims>& end_coords) {

  // The bounding box within which we will search for leaf entries
  ElemBoundingBox search_bbox(start_coords, end_coords);
    
  std::vector<std::reference_wrapper<const PayloadT>> found_entries;

  if(m_root_slot != NodePool::INVALID_SLOT) {
    search_overlapping_leaf_and_add(m_root_slot, search_bbox, found_entries);
  }
  
  return found_entries;
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::search_overlapping_leaf_and_add(std::size_t node_slot, const ElemBoundingBox& bbox,
									std::vector<std::reference_wrapper<const PayloadT>>& dest) {

  TreeNode& node = m_nodes[node_slot];

  if(node.level == 0) {

    // This is a leaf node; check if the bounding boxes of its data entries overlap ...
    for(std::size_t i = 0; i < node.num_children; i++) {      
      TreeEntry& entry = m_entries[node.child_slots[i]];

      // ... and, if they do, add them to `dest`
      if(bbox.overlaps(entry)) {
	dest.push_back(entry.payload);
      }
    }        
  }
  else {

    // This is an internal node in the tree; check if the bounding boxes of its children
    // overlap the search `bbox` ...
    for(std::size_t i = 0; i < node.num_children; i++) {
      std::size_t child_node_slot = node.child_slots[i];
      TreeNode& child_node = m_nodes[child_node_slot];

      // ... and, if they do, continue recursively
      if(bbox.overlaps(child_node)) {
	search_overlapping_leaf_and_add(child_node_slot, bbox, dest);
      }
    }
  }
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::GetAllEntries(std::vector<PayloadT>& dest) {  
  for(TreeEntry& cur_entry : m_entries) {
    dest.push_back(cur_entry.payload);
  }
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::DumpJSONTreeStructure(std::filesystem::path outpath) {
    
  std::ofstream outstream;
  outstream.open(outpath);  

  // Beginning of enclosing list of objects
  outstream << "[";

  // Dump everything
  if(m_root_slot != NodePool::INVALID_SLOT) {
    dump_JSON(m_root_slot, outstream, true);
  }

  // Close list with newline at the end of the file
  outstream << "\n]\n";
  
  outstream.close();
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::dump_JSON(std::size_t node_slot, std::ostream& stream, bool first_obj) {

  TreeNode& node = m_nodes[node_slot];
  
  if(node.level == 0) {

    // This is a leaf node; dump its own bounding box ...
    dump_node_JSON(node, stream, first_obj);
    
    // ... and all of its entries
    for(std::size_t i = 0; i < node.num_children; i++) {      
      TreeEntry& entry = m_entries[node.child_slots[i]];
      dump_entry_JSON(entry, stream, false);
    }
  }
  else {

    // This is an internal node in the tree; dump its own bounding box ...
    dump_node_JSON(node, stream, first_obj);

    // ... and continue with its children
    for(std::size_t i = 0; i < node.num_children; i++) {
      dump_JSON(node.child_slots[i], stream, false);
    }
  }  
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::dump_node_JSON(TreeNode& node, std::ostream& stream, bool first_obj) {

  if(first_obj) {
    stream << "\n{\n";  // begin object
  }
  else {
    stream << ",\n{\n"; // begin object
  }

  stream << "\"type\": \"Node\",\n";
  stream << "\"start_coords\": " << node.start_coords << ",\n";
  stream << "\"end_coords\": " << node.end_coords << ",\n";
  stream << "\"shape\": " << node.shape << ",\n";
  stream << "\"level\": " << node.level << "\n";

  // // --------------
  // stream << "\"num_children\": " << node.num_children << ",\n";
  // stream << "\"child_slots\": ";
  // for(std::size_t cur_child_slot: node.child_slots) {
  //   stream << cur_child_slot << ", ";
  // }  
  // stream << "\n";
  // // --------------
  
  stream << "}";  // end object
}

template <typename CoordT, std::size_t dims, class PayloadT>
void RStarTree<CoordT, dims, PayloadT>::dump_entry_JSON(TreeEntry& entry, std::ostream& stream, bool first_obj) {

  if(first_obj) {
    stream << "\n{\n";  // begin object
  }
  else {
    stream << ",\n{\n"; // begin object
  }

  stream << "\"type\": \"Entry\",\n";
  stream << "\"start_coords\": " << entry.start_coords << ",\n";
  stream << "\"end_coords\": " << entry.end_coords << ",\n";
  stream << "\"shape\": " << entry.shape << "\n";
  
  stream << "}";  // end object
}

// template <class IndexT, class PayloadT, std::size_t dims>
// void RStarTree<IndexT, PayloadT, dims>::Rebuild() {

//   // implement STL algorithm to rebuild tree

//   // 1) get a list of all TreeEntries from the memory pool
//   // 2) reset both memory pools
//   // 3) sort the list into groups such that each group goes into one TreeEntry
//   // 4) rebuild the new tree structure one group at a time:
//   //    4.1) add the elements in one group to subsequent entries in the entry pool
//   //    4.2) insert the corresponding nodes into the node pool
  
//   // -> also need to reorder the memory pools  
// }
