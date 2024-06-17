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
T& MemoryPool<T>::get_empty_slot() {

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
  return m_data[slot_ind];
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
template <class IndexT, typename PayloadIndT, typename NodeIndT, std::size_t MAX_NODESIZE>
Node<IndexT, PayloadIndT, NodeIndT, MAX_NODESIZE>::Node() : is_leaf(false), payload_ind(0), num_child_nodes(0) {
  child_inds.fill(0);
}

template <class IndexT, typename PayloadIndT, typename NodeIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, PayloadIndT, NodeIndT, MAX_NODESIZE>::mark_as_leaf(PayloadIndT payload_ind) {
  this -> payload_ind = payload_ind;
  is_leaf = true;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, typename PayloadIndT, typename NodeIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, PayloadIndT, NodeIndT, MAX_NODESIZE>::mark_as_internal() {
  this -> payload_ind = 0;
  is_leaf = false;
  num_child_nodes = 0;
  child_inds.fill(0);
}

template <class IndexT, typename PayloadIndT, typename NodeIndT, std::size_t MAX_NODESIZE>
void Node<IndexT, PayloadIndT, NodeIndT, MAX_NODESIZE>::add_child(NodeIndT child_ind) {
  assert(num_child_nodes < MAX_NODESIZE);     
  child_inds[num_child_nodes] = child_ind;
  num_child_nodes++;
}

// =======================

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
RTree<IndexT, PayloadT, dims, max_elems>::RTree(std::size_t init_slot_size) : m_nodes(init_slot_size), m_data(init_slot_size) { }

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
void RTree<IndexT, PayloadT, dims, max_elems>::AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind) {

}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
void RTree<IndexT, PayloadT, dims, max_elems>::Rebuild() {

}

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
// const PayloadT& RTree<IndexT, PayloadT, dims, max_elems>::Search(const IndexT& ind) {

// }

