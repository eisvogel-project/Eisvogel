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
  
  // fast cache entry lookup
  const PayloadT& get(const IndexT& elem);

  // inserts a new element into the cache, assuming there is an empty slot
  void insert_no_overwrite(const IndexT& index, const PayloadT& payload);

  // removes this element from the cache, creating a free slot
  const PayloadT& evict(const IndexT& elem);
  
  // assumes that the cache is full, evicts the oldest entry and returns it so that
  // it can be properly descoped
  const PayloadT& evict_oldest_full_cache();
  
  void print_old_to_new();
  void print_new_to_old();
  
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
