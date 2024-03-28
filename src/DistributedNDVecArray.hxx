#include "DistributedNDVecArray.hh"

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
CacheEntry<ArrayT, T, dims, vec_dims>& CacheEntry<ArrayT, T, dims, vec_dims>::operator=(std::size_t testval) {
  std::cout << "ultra special insertion assignment copy" << std::endl;
  return *this;
}

// --------------
