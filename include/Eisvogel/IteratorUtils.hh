#ifndef __ITERATOR_UTILS_HH
#define __ITERATOR_UTILS_HH

#include <cuchar>
#include "DenseNDArray.hh"

#include "Vector.hh"

// Iterator for static (i.e. known at compile time) number of dimensions
template <std::size_t vec_dims, class CallableT>
constexpr void loop_over_region(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
				CallableT&& worker) {
  Vector<std::size_t, vec_dims> cur;
  loop_over_dimensions<vec_dims - 1, vec_dims, CallableT>(begin, end, cur, worker);
}

template <std::size_t cur_dim, std::size_t vec_dims, class CallableT>
constexpr void loop_over_dimensions(const Vector<std::size_t, vec_dims>& begin, const Vector<std::size_t, vec_dims>& end,
				    Vector<std::size_t, vec_dims>& cur, CallableT&& worker) {

  static_assert((vec_dims > 0) && (cur_dim >= 0) && (cur_dim < vec_dims));

  for(cur[cur_dim] = begin[cur_dim]; cur[cur_dim] < end[cur_dim]; cur[cur_dim]++) {  
    if constexpr(cur_dim == 0) {
      worker(cur);
    }
    else {
      loop_over_dimensions<cur_dim - 1>(begin, end, cur, worker);
    }
  }  
}

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

#endif
