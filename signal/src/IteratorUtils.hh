#ifndef __ITERATOR_UTILS_HH
#define __ITERATOR_UTILS_HH

#include "NDArray.hh"

class IndexCounter {

private:

  IndexVector m_start, m_end, m_cur;
  std::size_t number_dims;

public:

  IndexCounter(const IndexVector& start, const IndexVector& end) : m_start(start), m_end(end), m_cur(start), 
								   number_dims(start.size()) { }
  IndexCounter& operator++() {
    for(size_t i = 0; i < number_dims; i++) {
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
    return m_cur(number_dims - 1) != m_end(number_dims - 1);
  }

  IndexVector::type& operator()(size_t ind) {return m_cur(ind);}

  auto begin() {return m_cur.begin();}
  auto end() {return m_cur.end();}

  IndexVector& index() {return m_cur;}
};

#endif
