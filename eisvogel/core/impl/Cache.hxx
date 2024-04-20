#include <cassert>

template <class IndexT, class PayloadT>
template <typename ... PayloadConstructorArgs>
Cache<IndexT, PayloadT>::Cache(std::size_t depth, PayloadConstructorArgs&& ... args) : m_depth(depth) {
  
  std::cout << "building cache with depth " << depth << std::endl;

  // allocate all required cache elements
  m_storage.reserve(depth); 
  for(std::size_t ind = 0; ind < depth; ind++) {
    m_storage.emplace_back(std::forward<PayloadConstructorArgs&&>(args)...);
  }

  // link them together
  CacheElement* prev_element = nullptr;
  for(CacheElement& cur_element : m_storage) {
    
    // link current entry to previous one ...
    cur_element.prev = prev_element;

    // ... and previous one to this one
    if(prev_element != nullptr) {
      prev_element -> next = &cur_element;
    }

    prev_element = &cur_element;
  }

  // keep track of front and back of the list
  m_oldest = &(m_storage.front());
  m_newest = &(m_storage.back());
}

template <class IndexT, class PayloadT>
template <class DataT>
void Cache<IndexT, PayloadT>::insert_no_overwrite(const IndexT& index, const DataT& payload) {

  PayloadT& insert_location = insert_ref_no_overwrite(index);
  insert_location = payload;
  
}

template <class IndexT, class PayloadT>
PayloadT& Cache<IndexT, PayloadT>::insert_ref_no_overwrite(const IndexT& index) {

  // an element with the new index must not already exist
  assert(!contains(index));
  
  // try to insert new element into the oldest slot
  CacheElement* insert_location = m_oldest;
  assert(!(insert_location -> occupied));

  // prepare the new payload ...
  insert_location -> index = index;
  insert_location -> occupied = true;

  // .. which now becomes the most-recently changed element in the cache
  mark_as_newest(insert_location);

  // update the index mapping
  m_accessor[index] = insert_location;
  
  return insert_location -> payload;
}

template <class IndexT, class PayloadT>
bool Cache<IndexT, PayloadT>::has_free_slot() {
  return m_accessor.size() < m_depth;
}

template <class IndexT, class PayloadT>
bool Cache<IndexT, PayloadT>::contains(const IndexT& elem) {
  return m_accessor.contains(elem);
}

template <class IndexT, class PayloadT>
std::vector<IndexT> Cache<IndexT, PayloadT>::contained_elements() {
  
  std::vector<IndexT> elements;
  elements.reserve(m_accessor.size());
  
  for (auto& [index, elem]: m_accessor) {
    elements.push_back(index);
  }
  return elements;
}

template <class IndexT, class PayloadT>
PayloadT& Cache<IndexT, PayloadT>::get(const IndexT& elem) {
    
  CacheElement* accessed_element = m_accessor.at(elem);

  // mark this as the most-recently accessed element
  mark_as_newest(accessed_element);
  
  return accessed_element -> payload;
}

template <class IndexT, class PayloadT>
PayloadT& Cache<IndexT, PayloadT>::evict_oldest_from_full_cache() {

  // this function assumes that the cache is full, i.e. that the oldest cache slot is actually occupied
  assert(!has_free_slot());
  
  // (don't need to move it, is already at the location of the oldest cache element)
  CacheElement* evicted_entry = m_oldest;
  evicted_entry -> occupied = false;

  // remove its entry from the index mapping
  m_accessor.erase(evicted_entry -> index);
  
  return evicted_entry -> payload;
}  

template <class IndexT, class PayloadT>
PayloadT& Cache<IndexT, PayloadT>::evict(const IndexT& elem) {

  // this function assumes that the element actually exists in the cache
  assert(contains(elem));
  
  CacheElement* evicted_entry = m_accessor.at(elem);
  evicted_entry -> occupied = false;

  // move it to the insert location at the "oldest" side of the cache
  mark_as_oldest(evicted_entry);

  // remove its entry from the index mapping
  m_accessor.erase(evicted_entry -> index);

  return evicted_entry -> payload;
}

// ------

template <class IndexT, class PayloadT>
void Cache<IndexT, PayloadT>::mark_as_oldest(CacheElement* element) {

  // this is already the oldest element, nothing needs to be done
  if(element -> prev == nullptr) {
    return;
  }

  // connect the elements to the left and to the right of this one
  element -> prev -> next = element -> next;
    
  if(element -> next != nullptr) {
    element -> next -> prev = element -> prev;
  }
  else {
    // `m_newest` points to `element`, make sure to update
    m_newest = element -> prev;
  }

  // move the element to the front
  element -> prev = nullptr;
  element -> next = m_oldest;
  m_oldest -> prev = element;

  // reseat pointer to the new element at the front of the queue
  m_oldest = element;
}

template <class IndexT, class PayloadT>
void Cache<IndexT, PayloadT>::mark_as_newest(CacheElement* element) {

  // this is already the newest element, nothing needs to be done
  if(element -> next == nullptr) {
    return;
  }

  // connect the elements to the left and to the right of this one
  element -> next -> prev = element -> prev;

  if(element -> prev != nullptr) {
    element -> prev -> next = element -> next;
  }
  else {    
    // `m_oldest` points to `element`, make sure to update
    m_oldest = element -> next;
  }

  // move the element to the back
  element -> next = nullptr;
  element -> prev = m_newest;
  m_newest -> next = element;

  // reseat pointer to the new element at the back of the queue 
  m_newest = element;  
}

template <class IndexT, class PayloadT>
template <class CallableT>
void Cache<IndexT, PayloadT>::loop_over_elements_old_to_new(CallableT&& worker) {

  CacheElement* cur_element = m_oldest;
  while(cur_element != nullptr) {
    
    if(cur_element -> occupied) {
      worker(cur_element -> payload);
    }
    cur_element = cur_element -> next;
  }  
}

template <class IndexT, class PayloadT>
template <class CallableT>
void Cache<IndexT, PayloadT>::loop_over_elements_new_to_old(CallableT&& worker) {

  CacheElement* cur_element = m_newest;
  while(cur_element != nullptr) {

    if(cur_element -> occupied) {
      worker(cur_element -> payload);
    }    
    cur_element = cur_element -> prev;
  }  
}
