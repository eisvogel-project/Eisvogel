#pragma once

#include <cuchar>
#include "Vector.hh"

// -----------------------------------------------

// Iterator for static (i.e. known at compile time) number of dimensions

namespace IteratorUtils {
    
  template <std::size_t cur_dim, std::size_t vec_dims, class NumT, class CallableT>
  constexpr void index_loop_over_elements_dimension(const Vector<NumT, vec_dims>& begin, const Vector<NumT, vec_dims>& end,
						    Vector<NumT, vec_dims>& cur, CallableT&& worker) {
    
    static_assert((vec_dims > 0) && (cur_dim >= 0) && (cur_dim < vec_dims));
    
    for(cur[cur_dim] = begin[cur_dim]; cur[cur_dim] < end[cur_dim]; cur[cur_dim]++) {  
      if constexpr(cur_dim == vec_dims - 1) {
	worker(cur);
      }
      else {
	// continue with iteration over the next-innermost index
	index_loop_over_elements_dimension<cur_dim + 1>(begin, end, cur, worker);
      }
    }  
  }

  template <std::size_t vec_dims, class NumT, class CallableT>
  constexpr void index_loop_over_elements(const Vector<NumT, vec_dims>& begin, const Vector<NumT, vec_dims>& end,
					  CallableT&& worker) requires(vec_dims == 2) {  
    Vector<NumT, vec_dims> cur;
    for(cur[0] = begin[0]; cur[0] < end[0]; cur[0]++) {
      for(cur[1] = begin[1]; cur[1] < end[1]; cur[1]++) {
	worker(cur);
      }
    }
  }
  
  template <std::size_t vec_dims, class NumT, class CallableT>
  constexpr void index_loop_over_elements(const Vector<NumT, vec_dims>& begin, const Vector<NumT, vec_dims>& end,
					  CallableT&& worker) {
    Vector<NumT, vec_dims> cur;
    index_loop_over_elements_dimension<0, vec_dims, NumT, CallableT>(begin, end, cur, worker);
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims, class CallableT>
  constexpr void index_loop_over_array_elements(const ArrayT<T, dims, vec_dims>& arr, CallableT&& worker) {
    Vector<std::size_t, dims> begin(0);
    Vector<std::size_t, dims> end = arr.GetShape();
    index_loop_over_elements(begin, end, worker);
  }
  
  // -----------------------------------------------
  
  // Iterator over chunks with certain size
  template <std::size_t cur_dim, std::size_t vec_dims, class NumT, class CallableT>
  constexpr void index_loop_over_chunks_dimension(const Vector<NumT, vec_dims>& begin, const Vector<NumT, vec_dims>& end,
						  const Vector<NumT, vec_dims>& chunk_size,
						  Vector<NumT, vec_dims>& chunk_begin, Vector<NumT, vec_dims>& chunk_end,
						  CallableT&& worker) {
    
    static_assert((vec_dims > 0) && (cur_dim >= 0) && (cur_dim < vec_dims));
    
    chunk_begin[cur_dim] = begin[cur_dim];
    while(chunk_begin[cur_dim] < end[cur_dim]) {
      chunk_end[cur_dim] = std::min(chunk_begin[cur_dim] + chunk_size[cur_dim], end[cur_dim]);
      
      if constexpr(cur_dim == vec_dims - 1) {
	worker(std::forward<Vector<NumT, vec_dims>>(chunk_begin), std::forward<Vector<NumT, vec_dims>>(chunk_end));
      }
      else {
	index_loop_over_chunks_dimension<cur_dim + 1>(begin, end, chunk_size, chunk_begin, chunk_end, worker);
      }
      
    chunk_begin[cur_dim] = chunk_end[cur_dim];
    }  
  }
    
  template <std::size_t vec_dims, class NumT, class CallableT>
  constexpr void index_loop_over_chunks(const Vector<NumT, vec_dims>& begin, const Vector<NumT, vec_dims>& end,
					const Vector<NumT, vec_dims>& chunk_size, CallableT&& worker) {
    Vector<NumT, vec_dims> chunk_begin;
    Vector<NumT, vec_dims> chunk_end;
    index_loop_over_chunks_dimension<0, vec_dims, NumT, CallableT>(begin, end, chunk_size, chunk_begin, chunk_end, worker);
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims, class CallableT>
  constexpr void index_loop_over_array_chunks(const ArrayT<T, dims, vec_dims>& arr,
					      const Vector<std::size_t, dims>& chunk_size, CallableT&& worker) {
    Vector<std::size_t, dims> begin(0u);
    Vector<std::size_t, dims> end = arr.GetShape();
    index_loop_over_chunks(begin, end, chunk_size, worker);
  }
}
