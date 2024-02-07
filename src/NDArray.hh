#ifndef __NDARRAY_HH
#define __NDARRAY_HH

#include <array>

// ======================================================
// General n-dimensional array
// ======================================================

template <class T, std::size_t dims>
class NDArray {

public:
  using shape_t = std::array<std::size_t, dims>;

public:
  NDArray() { }
  NDArray(const shape_t& shape) : m_shape(shape) { }
  virtual ~NDArray() = 0;

  const shape_t& shape() const {
    return m_shape;
  };

  std::size_t shape(std::size_t dim) const {
    return m_shape[dim];
  }

protected:
  shape_t m_shape = {};
};

template <class T, std::size_t dims>
inline NDArray<T, dims>::~NDArray() { }

#endif
