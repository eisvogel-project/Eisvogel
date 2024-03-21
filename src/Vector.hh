#ifndef __VECTOR_HH
#define __VECTOR_HH

#include <array>
#include <algorithm>
#include <stdexcept>
#include <execution>
#include <span>

// fixed-length vector
template <class T, std::size_t vec_dims>
class Vector {

private:
  using data_t = std::array<T, vec_dims>;

public:

  Vector() : m_data() { }

  Vector(std::span<T, vec_dims> view) {
    std::copy_n(std::execution::unseq, view.begin(), vec_dims, m_data.begin());
  }
  
  template <typename ... Values>
  Vector(Values ... values) requires(sizeof...(Values) == vec_dims) : m_data({values...}) { }
  
  T& operator[](std::size_t ind) {
    return m_data[ind];
  }

  const T& operator[](std::size_t ind) const {
    return m_data[ind];
  }  
  
  // includes bounds-checking
  T& at(std::size_t ind) {
    if(ind >= vec_dims) {
      throw std::out_of_range("");
    }
    return m_data[ind];
  }

  // arithmetic operations
  friend Vector<T, vec_dims> operator+(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::plus<>());
    return result;
  }

  Vector<T, vec_dims>& operator+=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::plus<>());
    return *this;
  }

  Vector<T, vec_dims>& operator*=(const T& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), m_data.begin(),
		   std::bind(std::multiplies<>(), std::placeholders::_1, rhs));
    return *this;
  }
  
  friend Vector<T, vec_dims> operator-(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::minus<>());
    return result;
  }

  Vector<T, vec_dims>& operator-=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::minus<>());
    return *this;
  }

  // iterators
  auto begin() {return m_data.begin();}
  auto cbegin() {return m_data.cbegin();}
  auto begin() const {return m_data.cbegin();}
  auto end() {return m_data.end();}
  auto cend() {return m_data.cend();}
  auto end() const {return m_data.cend();}
  
private:
  data_t m_data;

};

#endif
