#pragma once

#include <array>
#include <algorithm>
#include <stdexcept>
#include <execution>
#include <span>
#include <concepts>

#include "VectorView.hh"
#include "Serialization.hh"

// fixed-length vector
template <class T, std::size_t vec_dims> requires arithmetic<T>
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

  Vector(Vector<T, vec_dims>&& other) {
    std::copy_n(std::execution::unseq, other.m_data.begin(), vec_dims, m_data.begin());
  }
  
  // copy assignment
  Vector<T, vec_dims>& operator=(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.m_data.begin(), vec_dims, m_data.begin());
    return *this;
  }

  Vector<T, vec_dims>& operator=(const T& other) {
    std::fill_n(std::execution::unseq, m_data.begin(), vec_dims, other);
    return *this;
  }
  
  // Element access
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
  // type-casting
  // -----------------------------

  // TODO: make this private and collect the implementation in `Vector.hxx`
  template <typename DestT, typename SrcT, std::size_t... idxs>
  auto array_as_type(const std::array<SrcT, vec_dims>& src, std::index_sequence<idxs...>) const {
    return std::array<DestT, vec_dims>{{static_cast<DestT>(src[idxs])...}};
  }

  // returns a copy of this vector, with its values static_cast'ed to the destination type `DestT`
  template <class DestT>
  Vector<DestT, vec_dims> as_type() const {
    return Vector<DestT, vec_dims>(array_as_type<DestT>(m_data, std::make_index_sequence<vec_dims>()));
  }
  
  // -----------------------------
  // logic operations
  // -----------------------------

  friend bool operator==(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    for(std::size_t ind = 0; ind < vec_dims; ind++) {
      if(lhs[ind] != rhs[ind]) {
	return false;
      }
    }
    return true;
  }

  friend bool operator!=(const Vector<T, vec_dims>& lhs, const Vector<T, vec_dims>& rhs) {
    return !(lhs == rhs);
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

  // printing
  friend std::ostream& operator<<(std::ostream& stream, const Vector<T, vec_dims>& vec) {
    std::cout << "[ ";
    for(const T& cur: vec) {
      std::cout << cur << " ";
    }
    std::cout << "]";
    return stream;
  }
  
private:
  data_t m_data;

};

// Some useful utilities
namespace VectorUtils {

  template <typename T, std::size_t vec_dims>
  T inner(const Vector<T, vec_dims>& vec_a, const Vector<T, vec_dims>& vec_b);

  template <std::integral T1, std::integral T2, std::size_t vec_dims>
  Vector<T1, vec_dims> min(const Vector<T1, vec_dims>& vec_a, const Vector<T2, vec_dims>& vec_b);

  template <std::integral T1, std::integral T2, std::size_t vec_dims>
  Vector<T1, vec_dims> max(const Vector<T1, vec_dims>& vec_a, const Vector<T2, vec_dims>& vec_b);
  
}

// Some type shortcuts
template <typename T>
using Vector2D = Vector<T, 2>;

template <typename T>
using Vector3D = Vector<T, 3>;

template <typename T>
using Vector4D = Vector<T, 4>;

// Various kinds of vectors (for different semantics)
template <typename T>
struct RZVector : public Vector2D<T> {
  using Vector2D<T>::Vector2D;
  RZVector(Vector2D<T>&& other) : Vector2D<T>(std::forward<Vector2D<T>>(other)) { }
  RZVector(Vector2D<T>& other) : Vector2D<T>(std::forward<Vector2D<T>>(other)) { }

  const T& r() const { return this -> operator[](0); };
  const T& z() const { return this -> operator[](1); };

  T& r() { return this -> operator[](0); };
  T& z() { return this -> operator[](1); };

  static const T& r(const Vector2D<T>& vec) { return vec[0]; }
  static const T& z(const Vector2D<T>& vec) { return vec[1]; }

  static T& r(Vector2D<T>& vec) { return vec[0]; }
  static T& z(Vector2D<T>& vec) { return vec[1]; }
};

template <typename T>
struct ZRVector : public Vector2D<T> {
  using Vector2D<T>::Vector2D;
  ZRVector(Vector2D<T>&& other) : Vector2D<T>(std::forward<Vector2D<T>>(other)) { }
  ZRVector(Vector2D<T>& other) : Vector2D<T>(std::forward<Vector2D<T>>(other)) { }

  const T& z() const { return this -> operator[](0); };
  const T& r() const { return this -> operator[](1); };

  T& z() { return this -> operator[](0); };
  T& r() { return this -> operator[](1); };

  static const T& z(const Vector2D<T>& vec) { return vec[0]; }
  static const T& r(const Vector2D<T>& vec) { return vec[1]; }

  static T& z(Vector2D<T>& vec) { return vec[0]; }
  static T& r(Vector2D<T>& vec) { return vec[1]; }
};

template <typename T>
struct RZTVector : public Vector3D<T> {
  using Vector3D<T>::Vector3D;  // Inherit all constructors
  RZTVector(Vector3D<T>&& other) : Vector3D<T>(std::forward<Vector3D<T>>(other)) { }
  RZTVector(const Vector3D<T>& other) : Vector3D<T>(std::forward<const Vector3D<T>>(other)) { }
  RZTVector(const RZVector<T>& rz_vec, const T& t_val) : Vector3D<T>{rz_vec.r(), rz_vec.z(), t_val} { }

  const T& r() const { return this -> operator[](0); };
  const T& z() const { return this -> operator[](1); };
  const T& t() const { return this -> operator[](2); };

  T& r() { return this -> operator[](0); };
  T& z() { return this -> operator[](1); };
  T& t() { return this -> operator[](2); };

  static const T& r(const Vector3D<T>& vec) { return vec[0]; }
  static const T& z(const Vector3D<T>& vec) { return vec[1]; }
  static const T& t(const Vector3D<T>& vec) { return vec[2]; }

  RZVectorView<T> rz_view() { return RZVectorView<T>( this -> begin()); };
};

template <typename T>
struct TZRVector : public Vector3D<T> {
  using Vector3D<T>::Vector3D;
  TZRVector(Vector3D<T>&& other) : Vector3D<T>(std::forward<Vector3D<T>>(other)) { }

  const T& t() const { return this -> operator[](0); };
  const T& z() const { return this -> operator[](1); };
  const T& r() const { return this -> operator[](2); };

  T& t() { return this -> operator[](0); };
  T& z() { return this -> operator[](1); };
  T& r() { return this -> operator[](2); };
};

template <typename T>
struct XYZTVector : public Vector4D<T> {
  using Vector4D<T>::Vector4D;
  XYZTVector(Vector4D<T>&& other) : Vector4D<T>(std::forward<Vector4D<T>>(other)) { }

  const T& x() const { return this -> operator[](0); };
  const T& y() const { return this -> operator[](1); };
  const T& z() const { return this -> operator[](2); };
  const T& t() const { return this -> operator[](3); };

  T& x() { return this -> operator[](0); };
  T& y() { return this -> operator[](1); };
  T& z() { return this -> operator[](2); };
  T& t() { return this -> operator[](3); };
};

template <typename T>
struct XYZVector : Vector3D<T> {
  using Vector3D<T>::Vector3D;
  XYZVector(Vector3D<T>&& other) : Vector3D<T>(std::forward<Vector3D<T>>(other)) { }

  const T& x() const { return this -> operator[](0); };
  const T& y() const { return this -> operator[](1); };
  const T& z() const { return this -> operator[](2); };

  T& x() { return this -> operator[](0); };
  T& y() { return this -> operator[](1); };
  T& z() { return this -> operator[](2); };
};

struct XYZTCoordVector : XYZTVector<scalar_t> {
  using XYZTVector<scalar_t>::XYZTVector;
};
struct XYZTFieldVector : XYZTVector<scalar_t> {
  using XYZTVector<scalar_t>::XYZTVector;
};
struct XYZTIndexVector : XYZTVector<std::size_t> {
  using XYZTVector<std::size_t>::XYZTVector;
};

struct XYZCoordVector : XYZVector<scalar_t> {
  using XYZVector<scalar_t>::XYZVector;
};
struct XYZFieldVector : XYZVector<scalar_t> {
  using XYZVector<scalar_t>::XYZVector;
};

struct RZCoordVector : RZVector<scalar_t> {
  using RZVector<scalar_t>::RZVector;
};
struct RZFieldVector : RZVector<scalar_t> {
  using RZVector<scalar_t>::RZVector;
};
struct RZIndexVector : RZVector<std::size_t> {
  using RZVector<std::size_t>::RZVector;
};

struct ZRCoordVector : ZRVector<scalar_t> {
  using ZRVector<scalar_t>::ZRVector;
};
struct ZRFieldVector : ZRVector<scalar_t> {
  using ZRVector<scalar_t>::ZRVector;
};
struct ZRIndexVector : ZRVector<std::size_t> {
  using ZRVector<std::size_t>::ZRVector;
};

struct RZTCoordVector : RZTVector<scalar_t> {
  using RZTVector<scalar_t>::RZTVector;
};
struct RZTIndexVector : RZTVector<std::size_t> {
  using RZTVector<std::size_t>::RZTVector;
};
struct RZTSignedIndexVector : RZTVector<int> {
  using RZTVector<int>::RZTVector;
};

struct TZRCoordVector : TZRVector<scalar_t> {
  using TZRVector<scalar_t>::TZRVector;
};
struct TZRIndexVector : TZRVector<std::size_t> {
  using TZRVector<std::size_t>::TZRVector;
};

#include "Vector.hxx"
