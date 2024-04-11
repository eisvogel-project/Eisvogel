#pragma once

#include <span>
#include <iterator>
#include <algorithm>
#include "Eisvogel/Common.hh"

// forward declaration
template<typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

template <class T, std::size_t vec_dims> requires arithmetic<T>
class Vector;

template <typename T, std::size_t vec_dims>
struct VectorView : public std::span<T, vec_dims> {

  // create a view to a memory location described by the iterator `It`
  template <std::forward_iterator It>
  constexpr VectorView(It first) : std::span<T, vec_dims>(first, vec_dims) { }
  
  // create a view to a given vector
  constexpr VectorView(Vector<T, vec_dims>& vec) : VectorView(vec.begin()) { }
  
  // copy constructor
  constexpr VectorView(const VectorView& other) : std::span<T, vec_dims>(other) { }
  
  VectorView& operator=(const Vector<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.cbegin(), vec_dims, this -> begin());
    return *this;
  }

  VectorView& operator=(const VectorView<T, vec_dims>& other) {
    std::copy_n(std::execution::unseq, other.begin(), vec_dims, this -> begin());
    return *this;
  }
};

// Some type shortcuts
template <typename T>
using Vector2DView = VectorView<T, 2>;

template <typename T>
using Vector3DView = VectorView<T, 3>;

template <typename T>
using Vector4DView = VectorView<T, 4>;

// Various kinds of vector views (for different semantics)
template <typename T>
struct RZVectorView : public Vector2DView<T> {
  using Vector2DView<T>::Vector2DView;
  RZVectorView(const Vector2DView<T>& other) : Vector2DView<T>(other) { }

  static const T& r(const RZVectorView& vec) { return vec[0]; }
  static const T& z(const RZVectorView& vec) { return vec[1]; }

  static T& r(RZVectorView& vec) { return vec[0]; }
  static T& z(RZVectorView& vec) { return vec[1]; }
};

// More typedefs (can later turn them into their own types if needed)
using RZCoordVectorView = RZVectorView<scalar_t>;
using RZFieldVectorView = RZVectorView<scalar_t>;
