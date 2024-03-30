#pragma once

#include <vector>
#include <unordered_map>

template <class IndexT, class PayloadT>
class Cache {
  
public:  
  
  template <typename ... PayloadArgs>
  Cache(std::size_t depth, PayloadArgs&& ... args);

  bool has_free_slot();

  bool contains(const IndexT& elem);
  std::vector<IndexT> contained_elements();
  
  // fast cache entry lookup
  PayloadT& get(const IndexT& elem);

  // inserts a new element into the cache, assuming there is an empty slot
  template <class DataT>
  void insert_no_overwrite(const IndexT& index, const DataT& payload);

  // reserves a new element in the cache and returns a reference to that clients can directly write to it
  PayloadT& insert_ref_no_overwrite(const IndexT& index);

  // removes this element from the cache, creating a free slot and returning a reference to the element at this location
  // cache slot can be overwritten at the next insert call
  // is now the responsibility of the caller to do anything that needs to be done to not lose the information
  PayloadT& evict(const IndexT& elem);
  
  // assumes that the cache is full, evicts the oldest entry and returns it so that
  // it can be properly descoped
  // cache slot can be overwritten at the next insert call
  // is now the responsibility of the caller to do anything that needs to be done to not lose the information
  PayloadT& evict_oldest_from_full_cache();
  
  template <class CallableT>
  void loop_over_elements_old_to_new(CallableT&& worker);

  template <class CallableT>
  void loop_over_elements_new_to_old(CallableT&& worker);
  
private:

  struct CacheElement {
    
    template <typename ... PayloadConstructorArgs>
    CacheElement(PayloadConstructorArgs&& ... args) :
      payload(std::forward<PayloadConstructorArgs&&>(args)...), occupied(false), next(nullptr), prev(nullptr) { };
    
    PayloadT payload;
    IndexT index;
    bool occupied;
    
    CacheElement* next;
    CacheElement* prev;
  }; 

  // moves to the "oldestmost" position of the list
  void mark_as_oldest(CacheElement* element);

  // moves to the "newestmost" position of the list
  void mark_as_newest(CacheElement* element);

private:

  std::size_t m_depth;
  
  std::unordered_map<IndexT, CacheElement*> m_accessor;
  std::vector<CacheElement> m_storage;
  CacheElement* m_oldest;
  CacheElement* m_newest;
};

#include "Cache.hxx"
