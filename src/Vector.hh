#ifndef __VECTOR_HH
#define __VECTOR_HH

#include <array>
#include <algorithm>
#include <stdexcept>
#include <execution>
#include <span>

#include "Serialization.hh"

// fixed-length vector
template <class T, std::size_t vec_dims>
class Vector {

private:
  using data_t = std::array<T, vec_dims>;

private:
  friend struct stor::Traits<Vector<T, vec_dims>>;
  
public:

  Vector() : m_data() { }

  Vector(std::span<T, vec_dims> view) {
    std::copy_n(std::execution::unseq, view.begin(), vec_dims, m_data.begin());
  }
  
  template <typename ... Values>
  Vector(Values ... values) requires(sizeof...(Values) == vec_dims) : m_data({values...}) { }

  Vector(const T& value) {
    std::fill_n(std::execution::unseq, m_data.begin(), vec_dims, value);
  }

  Vector(const data_t&& data) {
    std::copy_n(std::execution::unseq, data.begin(), vec_dims, m_data.begin());
  }

  // copy constructor
  Vector(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.m_data.begin(), vec_dims, m_data.begin());
  }
  
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

  // -----------------------------
  // arithmetic operations
  // -----------------------------

  // element-wise operations
  friend Vector<T, vec_dims> operator+(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::plus<>());
    return result;
  }

  friend Vector<T, vec_dims> operator-(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::minus<>());
    return result;
  }
  
  friend Vector<T, vec_dims> operator*(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::multiplies<>());
    return result;
  }

  friend Vector<T, vec_dims> operator/(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(),
		   std::divides<>());
    return result;
  }  

  // element-wise operations (in-place)
  Vector<T, vec_dims>& operator+=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::plus<>());
    return *this;
  }

  Vector<T, vec_dims>& operator-=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::minus<>());
    return *this;
  }

  Vector<T, vec_dims>& operator*=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::multiplies<>());
    return *this;
  }

  Vector<T, vec_dims>& operator/=(const Vector<T, vec_dims>& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(),
		   std::divides<>());
    return *this;
  }

  // scalar operations
  friend Vector<T, vec_dims> operator+(const Vector<T, vec_dims>& lhs, const T& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), result.m_data.begin(),
		   std::bind(std::plus<>(), std::placeholders::_1, rhs));
    return result;
  }  

  friend Vector<T, vec_dims> operator-(const Vector<T, vec_dims>& lhs, const T& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), result.m_data.begin(),
		   std::bind(std::minus<>(), std::placeholders::_1, rhs));
    return result;
  }  

  friend Vector<T, vec_dims> operator*(const Vector<T, vec_dims>& lhs, const T& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), result.m_data.begin(),
		   std::bind(std::multiplies<>(), std::placeholders::_1, rhs));
    return result;
  }  

  friend Vector<T, vec_dims> operator/(const Vector<T, vec_dims>& lhs, const T& rhs) {
    Vector<T, vec_dims> result;
    std::transform(std::execution::unseq, lhs.m_data.begin(), lhs.m_data.end(), result.m_data.begin(),
		   std::bind(std::divides<>(), std::placeholders::_1, rhs));
    return result;
  }  
  
  // scalar operations (in-place)
  Vector<T, vec_dims> operator+=(const T& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), m_data.begin(),
		   std::bind(std::plus<>(), std::placeholders::_1, rhs));
    return *this;
  }

  Vector<T, vec_dims> operator-=(const T& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), m_data.begin(),
		   std::bind(std::minus<>(), std::placeholders::_1, rhs));
    return *this;
  }
  
  Vector<T, vec_dims>& operator*=(const T& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), m_data.begin(),
		   std::bind(std::multiplies<>(), std::placeholders::_1, rhs));
    return *this;
  }

  Vector<T, vec_dims>& operator/=(const T& rhs) {
    std::transform(std::execution::unseq, m_data.begin(), m_data.end(), m_data.begin(),
		   std::bind(std::divides<>(), std::placeholders::_1, rhs));
    return *this;
  }
  
  // iterators
  auto begin() {return m_data.begin();}
  auto cbegin() const {return m_data.cbegin();}
  auto begin() const {return m_data.cbegin();}
  auto end() {return m_data.end();}
  auto cend() const {return m_data.cend();}
  auto end() const {return m_data.cend();}
  
private:
  data_t m_data;

};

namespace stor {

  template <typename T, std::size_t vec_dims>
  struct Traits<Vector<T, vec_dims>> {
    using type = Vector<T, vec_dims>;
    using data_t = typename type::data_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<data_t>::serialize(stream, val.m_data);
    }

    static type deserialize(std::iostream& stream) {
      data_t data = Traits<data_t>::deserialize(stream);
      return Vector<T, vec_dims>(std::move(data));
    }
  };    
}

// Some type shortcuts
template <typename T>
using Vector2D = Vector<T, 2>;

template <typename T>
using Vector3D = Vector<T, 3>;

template <typename T>
using Vector4D = Vector<T, 4>;

using IndexVector2D = Vector2D<std::size_t>;
using IndexVector3D = Vector3D<std::size_t>;
using IndexVector4D = Vector4D<std::size_t>;

#endif
