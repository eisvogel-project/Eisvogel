#include <algorithm>

template <class T>
MemoryPool<T>::MemoryPool(std::size_t init_size) : m_data(init_size), m_slots_free(init_size) {

  // Mark all slots as free
  std::iota(m_slots_free.begin(), m_slots_free.end(), 0);
}

template <class T>
T& MemoryPool<T>::operator[](IndT ind) {
  return m_data[ind];
}

template <class T>
const T& MemoryPool<T>::operator[](IndT ind) const {
  return m_data[ind];
}

template <class T>
MemoryPool<T>::IndT MemoryPool<T>::get_empty_slot() {

  // Grow the pool if needed
  if(m_slots_free.size() == 0) {
    grow();
  }

  // Get the index of the next free slot
  // ----------------------------------------------------------------------
  // Note: at the moment, this fills the memory slots from the back, which
  // does not play optimally with cache lines
  // ----------------------------------------------------------------------
  // TODO: turn this into a doubly-linked list where removal of elements
  // from both ends is cheap
  // ----------------------------------------------------------------------
  IndT slot_ind = m_slots_free.back();
  m_slots_free.pop_back();
  return slot_ind;
}

template <class T>
void MemoryPool<T>::grow() {

  // Double the size of the currently-allocated memory
  std::size_t current_size = m_data.size();
  std::size_t new_size = 2 * current_size;
  m_data.resize(new_size);

  // Mark the newly-created slots as free
  std::vector<IndT> new_slots_free(new_size - current_size);
  std::iota(new_slots_free.begin(), new_slots_free.end(), current_size);
  m_slots_free.insert(m_slots_free.end(), new_slots_free.begin(), new_slots_free.end());
}

template <class T>
bool MemoryPool<T>::is_free(IndT ind) {
  return std::find(m_slots_free.begin(), m_slots_free.end(), ind) != m_slots_free.end();
}

template <class T>
bool MemoryPool<T>::is_allocated(IndT ind) {
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

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::RTree(std::size_t init_slot_size) : m_nodes(init_slot_size), m_data(init_slot_size) {

  // Start with an empty root node
  m_root_node_ind = m_nodes.get_empty_slot();
  m_nodes[m_root_node_ind].mark_as_empty_leaf();
}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
void RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind) {

  // Take ownership of the `elem` and store it
  PayloadIndT payload_ind = m_data.get_empty_slot();
  m_data[payload_ind] = elem;

  // Build a tree node that points to it ...
  NodeIndT node_ind = m_nodes.get_empty_slot();
  m_nodes[node_ind].mark_as_empty_leaf();

  // ... and make this the child of an existing node in the tree
  
}

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
// void RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::Rebuild() {

// }

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t MAX_NODESIZE>
// const PayloadT& RTree<IndexT, PayloadT, dims, MAX_NODESIZE>::Search(const IndexT& ind) {

// }

