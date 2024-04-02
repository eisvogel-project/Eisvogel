#pragma once

#include <cuchar>
#include "DenseNDArray.hh"
#include "Vector.hh"

// -----------------------------------------------

// Iterator for static (i.e. known at compile time) number of dimensions

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_array_elements(const ArrayT<T, dims, vec_dims>& arr, CallableT&& worker) {
  Vector<std::size_t, dims> begin(0);
  Vector<std::size_t, dims> end = arr.GetShape();
  index_loop_over_elements(begin, end, worker);
}

template <std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_elements(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
					CallableT&& worker) requires(vec_dims == 2) {  
  Vector<std::size_t, vec_dims> cur;
  for(cur[0] = begin[0]; cur[0] < end[0]; cur[0]++) {
    for(cur[1] = begin[1]; cur[1] < end[1]; cur[1]++) {
      worker(cur);
    }
  }
}

template <std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_elements(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
					CallableT&& worker) {
  Vector<std::size_t, vec_dims> cur;
  index_loop_over_elements_dimension<0, vec_dims, CallableT>(begin, end, cur, worker);
}

template <std::size_t cur_dim, std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_elements_dimension(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
						  Vector<std::size_t, vec_dims>& cur, CallableT&& worker) {

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

// -----------------------------------------------

// Iterator over chunks with certain size

template <template<typename, std::size_t, std::size_t> class ArrayT,
	  typename T, std::size_t dims, std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_array_chunks(const ArrayT<T, dims, vec_dims>& arr,
					    const Vector<std::size_t, dims>& chunk_size, CallableT&& worker) {
  Vector<std::size_t, dims> begin(0);
  Vector<std::size_t, dims> end = arr.GetShape();
  index_loop_over_chunks(begin, end, chunk_size, worker);
}

template <std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_chunks(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
				      const Vector<std::size_t, vec_dims>& chunk_size, CallableT&& worker) {
  Vector<std::size_t, vec_dims> chunk_begin;
  Vector<std::size_t, vec_dims> chunk_end;
  index_loop_over_chunks_dimension<0, vec_dims, CallableT>(begin, end, chunk_size, chunk_begin, chunk_end, worker);
}

template <std::size_t cur_dim, std::size_t vec_dims, class CallableT>
constexpr void index_loop_over_chunks_dimension(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
						const Vector<std::size_t, vec_dims>& chunk_size,
						Vector<std::size_t, vec_dims>& chunk_begin, Vector<std::size_t, vec_dims>& chunk_end,
						CallableT&& worker) {
  
  static_assert((vec_dims > 0) && (cur_dim >= 0) && (cur_dim < vec_dims));

  chunk_begin[cur_dim] = begin[cur_dim];
  while(chunk_begin[cur_dim] < end[cur_dim]) {
    chunk_end[cur_dim] = std::min(chunk_begin[cur_dim] + chunk_size[cur_dim], end[cur_dim]);

    if constexpr(cur_dim == vec_dims - 1) {
      worker(chunk_begin, chunk_end);
    }
    else {
      index_loop_over_chunks_dimension<cur_dim + 1>(begin, end, chunk_size, chunk_begin, chunk_end, worker);
    }

    chunk_begin[cur_dim] = chunk_end[cur_dim];
  }  
}


// -----------------------------------------------

// Iterator for dynamic (i.e. known at run time) number of dimensions
template <typename VecT> class VectorCounter {

private:

  VecT m_start, m_end, m_cur;
  std::size_t number_dims;

public:

  VectorCounter(const VecT& start, const VecT& end) : m_start(start), m_end(end), m_cur(start), 
						      number_dims(start.size()) { }
  VectorCounter& operator++() {
    for(std::size_t i = 0; i < number_dims; i++) {
      if(++m_cur(i) == m_end(i)) {
	if(i < number_dims - 1) {
	  m_cur(i) = m_start(i);
	}
      }
      else
	break;
    }

    return *this;
  }

  bool running() {
    return m_cur(number_dims - 1) < m_end(number_dims - 1);
  }

  VecT::type& operator()(std::size_t ind) {return m_cur(ind);}

  auto begin() {return m_cur.begin();}
  auto end() {return m_cur.end();}

  inline VecT& index() {return m_cur;}
};

using IndexCounter = VectorCounter<IndexVector>;
using GridCounter = VectorCounter<GridVector>;

// some additional utility functions
bool isInIndexRange(const IndexVector& inds, const IndexVector& start_inds, const IndexVector& stop_inds);
