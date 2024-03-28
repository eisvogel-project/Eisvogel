#pragma once

#include <unordered_map>

template <class IndexT, class PayloadT>
struct CacheEntry {

  template <typename ... ConstructorArgs>
  CacheEntry(ConstructorArgs&& ... args) : payload(args ...) { };
  
  PayloadT payload;
  IndexT index;

  CacheEntry<IndexT, PayloadT>* newer;
  CacheEntry<IndexT, PayloadT>* older;
};

template <class IndexT, class PayloadT>
class Cache {

  using cache_entry_t = CacheEntry<IndexT, PayloadT>;
  
public:  
  
  template <typename ... PayloadArgs>
  Cache(std::size_t depth, PayloadArgs&& ... args) {
    for(std::size_t ind = 0; ind < depth; ind++) {
      // build vector of CacheEntries and link them
    }
  }

  // fast cache entry lookup
  const PayloadT& get(IndexT& elem) {
    // cache_entry_t* accessed_element = accessor[elem]
    // move to front of list
    // return accessed_element -> payload;
  }

  bool has_free_slot() {
    // return accessor.size() < m_depth;
    return false;
  }

  void insert(const PayloadT& payload) {
    // cache_entry_t* insert_location = oldest;
    // insert_location -> payload = payload; // need to properly overload assignment to copy the things that need to be copied
    
    // insert into list
    // mark_as_newest(insert_location);
  }

  // removes oldest entry from the cache and returns it so that it can be properly descoped
  const PayloadT& remove_oldest() {
    // cache_entry_t* removed_entry = oldest;
    // return removed_entry -> payload;
  }
  
  // removes this element from the cache, creating a free slot
  void free(IndexT& elem) {
    // cache_entry_t* this_entry = accessor[elem];
    // accessor.erase(this_entry -> index); // remove from accessor
    // mark_as_oldest(this_entry);
  }

private:

  // moves to the "oldestmost" position of the list
  void mark_as_oldest(cache_entry_t* element);

  // moves to the "newestmost" position of the list
  void mark_as_newest(cache_entry_t* element);

private:

  // Cache size
  std::size_t m_depth;
  
  std::unordered_map<IndexT, cache_entry_t*> accessor;
  std::vector<cache_entry_t> storage;
  cache_entry_t* oldest;
  cache_entry_t* newest;
};
