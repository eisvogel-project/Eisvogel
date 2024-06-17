#include <algorithm>

template <class T, std::unsigned_integral IndT>
MemoryPool<T, IndT>::MemoryPool(std::size_t init_size) : m_init_size(init_size),
							 m_free_start(Slot::INVALID), m_free_end(Slot::INVALID),
							 m_alloc_start(Slot::INVALID), m_alloc_end(Slot::INVALID),
							 m_data(init_size), m_slots(init_size) {
  assert(m_init_size > 0);
  reset_slot_lists();
}

template <class T, std::unsigned_integral IndT>
void MemoryPool<T, IndT>::reset() {

  m_data.clear();
  m_data.resize(m_init_size);

  m_slots.clear();
  m_slots.resize(m_init_size);
  
  reset_slot_lists();
}

template <class T, std::unsigned_integral IndT>
void MemoryPool<T, IndT>::reset_slot_lists() {

  assert(m_data.size() == m_slots.size());

  std::size_t number_slots = m_slots.size();
  
  // Set the two end points of the doubly-linked lists:
  // Nothing is allocated yet ...
  m_alloc_start = Slot::INVALID; 
  m_alloc_end = Slot::INVALID;  

  // ... and everything is free
  m_free_start = 0;
  m_free_end = number_slots - 1;

  // Link up the list
  for(IndT i = 0; i < number_slots; i++) {
    m_slots[i].data_ind = i;
    m_slots[i].next_data_ind = i + 1;
    m_slots[i].prev_data_ind = i - 1;
  }
  m_slots.front().prev_data_ind = Slot::INVALID;
  m_slots.back().next_data_ind = Slot::INVALID;   
}

template <class T, std::unsigned_integral IndT>
T& MemoryPool<T, IndT>::operator[](IndT ind) {
  return m_data[ind];
}

template <class T, std::unsigned_integral IndT>
const T& MemoryPool<T, IndT>::operator[](IndT ind) const {
  return m_data[ind];
}

template <class T, std::unsigned_integral IndT>
IndT MemoryPool<T, IndT>::get_empty_slot_front() {

  // Grow the pool if needed
  if(m_free_start == Slot::Invalid) {
    grow();
  }

  // Get the index of the next free slot
  IndT alloc_slot_ind = m_free_start;

  // Update the lists
  
  if(m_alloc_start == Slot::INVALID) {
    assert(m_alloc_end == Slot::INVALID);
    m_alloc_start = alloc_slot_ind;
    m_alloc_end = alloc_slot_ind;

    m_slots[alloc_slot_ind].next_data_ind = Slot::INVALID;
    m_slots[alloc_slot_ind].prev_data_ind = Slot::INVALID;
  }
  else {  
    m_slots[m_alloc_end].next_data_ind = alloc_slot_ind;
    m_slots[alloc_slot_ind].prev_data_ind = m_alloc_end;
    m_alloc_end = alloc_slot_ind;
  }
  
  return alloc_slot_ind;
}

template <class T, std::unsigned_integral IndT>
void MemoryPool<T, IndT>::grow() {

  // Double the size of the currently-allocated memory
  std::size_t current_size = m_data.size();
  std::size_t new_size = 2 * current_size;
  m_data.resize(new_size);

  // Extend the slot list
  std::vector<Slot> new_slots_free(new_size - current_size);
  m_slots.insert(m_slots.end(), new_slots_free.begin(), new_slots_free.end());  
  for(IndT i = current_size; i < new_size; i++) {
    m_slots[i].data_ind = i;
    m_slots[i].next_data_ind = i + 1;
    m_slots[i].prev_data_ind = i - 1;
  }
  m_slots[current_size].prev_data_ind = m_free_end; // Points back to the element that was previously the last free slot
                                                    // (irrespective of whether it was Slot::INVALID or not)
  m_slots.back().next_data_ind = Slot::INVALID;

  if(m_free_end != Slot::INVALID) {
    m_slots[m_free_end].next_data_ind = current_size;
  } 

  if(m_free_start == Slot::INVALID) {
    m_free_start = current_size;
  }
  m_free_end = new_size - 1;
}

template <class T, std::unsigned_integral IndT>
void MemoryPool<T, IndT>::print_free() {

  IndT cur_free = m_free_start;
  
  while(true) {
    if(cur_free == Slot::INVALID) {
      break;
    }

    std::cout << cur_free << std::endl;
    cur_free = m_slots[cur_free].next_data_ind;
  }
}

template <class T, std::unsigned_integral IndT>
bool MemoryPool<T, IndT>::is_free(IndT ind) {
  return false; // TODO
}

template <class T, std::unsigned_integral IndT>
bool MemoryPool<T, IndT>::is_allocated(IndT ind) {
  return !is_free(ind);
}

// =======================

// Default constructor: mark everything as invalid
template <class IndexT, typename IndT, std::size_t MAX_NODESIZE>
Node<IndexT, IndT, MAX_NODESIZE>::Node() : is_leaf(true), num_child_nodes(0) {
  child_inds.fill(0);
}

template <class IndexT, typename IndT, std::size_t MAX_NODESIZE>
void Node<IndexT, IndT, MAX_NODESIZE>::mark_as_empty_leaf() {
  is_leaf = true;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, typename IndT, std::size_t MAX_NODESIZE>
void Node<IndexT, IndT, MAX_NODESIZE>::mark_as_empty_internal() {
  is_leaf = false;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, typename IndT, std::size_t MAX_NODESIZE>
void Node<IndexT, IndT, MAX_NODESIZE>::add_child(IndT child_ind) {
  assert(num_child_nodes < MAX_NODESIZE);     
  child_inds[num_child_nodes] = child_ind;
  num_child_nodes++;
}

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::RTree(std::size_t init_slot_size) : m_nodes(init_slot_size), m_entries(init_slot_size) {

  // Start with an empty root node
  m_root_node_ind = m_nodes.get_empty_slot();
  m_nodes[m_root_node_ind].mark_as_empty_leaf();
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::IndT RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::choose_subtree(IndT start_node) {
  
  return start_node;
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE, std::size_t MIN_NODESIZE>
void RTree<IndexT, PayloadT, dims, MAX_NODESIZE, MIN_NODESIZE>::AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind) {

  // Take ownership of the `elem` and store it
  IndT payload_ind = m_entries.get_empty_slot();
  m_entries[payload_ind].payload = elem;

  // Build a tree node that points to it ...
  IndT node_ind = m_nodes.get_empty_slot();
  m_nodes[node_ind].mark_as_empty_leaf();

  // ... and make this the child of an existing node in the tree
  
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

