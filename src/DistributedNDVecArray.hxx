#include "DistributedNDVecArray.hh"

// --------------

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims>
CacheEntry<ArrayT, T, dims, vec_dims>& CacheEntry<ArrayT, T, dims, vec_dims>::operator=(std::tuple<metadata_t&, chunk_t&, CacheStatus&> other) {

  std::cout << "ultra special insertion assignment copy" << std::endl;

  std::cout << "metadata copy" << std::endl;
  chunk_meta = std::get<metadata_t&>(other);

  std::cout << "NDVec array copy" << std::endl;
  chunk_data = std::get<chunk_t&>(other);
  
  op_to_perform = std::get<CacheStatus&>(other); 
  
  return *this;
}

// --------------
