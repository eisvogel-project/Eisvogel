#ifndef __NDARRAY_HH
#define __NDARRAY_HH

#include <vector>
#include <array>
#include <numeric>
#include <algorithm>

template <class T, std::size_t dims>
class NDArray {

protected:
  using shape_t = std::array<std::size_t, dims>;

public:
  NDArray(const shape_t& shape) : m_shape(shape) { }
  virtual ~NDArray() = 0;

  const shape_t& shape() const {
    return m_shape;
  };

  std::size_t shape(std::size_t dim) const {
    return m_shape[dim];
  }

protected:
  const shape_t m_shape = {};
};

template <class T, std::size_t dims>
inline NDArray<T, dims>::~NDArray() { }

// ======

template <class T, std::size_t dims>
class DenseNDArray : public NDArray<T, dims> {

public:
  using shape_t = typename NDArray<T, dims>::shape_t;

private:
  using stride_t = std::array<std::size_t, dims + 1>;
  using data_t = std::vector<T>;

public:

  using type = T;

  DenseNDArray(std::size_t size, const T& value) requires(dims == 1) : DenseNDArray(shape_t({size}), value) { }
  
  DenseNDArray(const shape_t& shape, const T& value) : NDArray<T, dims>(shape) {
    m_strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), m_strides.begin() + 1, std::multiplies<std::size_t>());
    m_data.resize(m_strides.back(), value);
  }

  DenseNDArray(const shape_t& shape, std::vector<T>&& data) : NDArray<T, dims>(shape) {
    m_strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), m_strides.begin() + 1, std::multiplies<std::size_t>());

    if(m_strides.back() != data.size())
      throw;

    m_data = std::move(data);
  }

  DenseNDArray(std::initializer_list<T>&& data) requires(dims == 1) :
    NDArray<T, 1>({data.size()}), m_strides({1, data.size()}), m_data(data.begin(), data.end()) { }

  DenseNDArray(std::vector<T>&& data) requires(dims == 1) : 
    NDArray<T, 1>({data.size()}), m_strides({1, data.size()}), m_data(data) { }

  DenseNDArray(const DenseNDArray<T, dims>& other) : 
    NDArray<T, dims>(other.m_shape), m_strides(other.m_strides), m_data(other.m_data) { }

  // Indexing with explicit pack of coordinates
  template <typename... Inds>
  const T& operator()(Inds... inds) const requires(sizeof...(Inds) == dims) {
    std::size_t flat_ind = 0, dim = 0;
    (..., (flat_ind += inds * m_strides[dim++]));
    return m_data[flat_ind];
  }

  template <typename... Inds>
  T& operator()(Inds... inds) requires(sizeof...(Inds) == dims) {
    std::size_t flat_ind = 0, dim = 0;
    (..., (flat_ind += inds * m_strides[dim++]));
    return m_data[flat_ind];
  }

  // Indexing with a vector that holds the coordinates
  const T& operator()(DenseNDArray<std::size_t, 1>& inds) const {
    if(inds.size() != dims)
      throw;
    
    std::size_t flat_ind = std::inner_product(inds.cbegin(), inds.cend(), m_strides.begin(), 0);
    return m_data[flat_ind];
  }

  T& operator()(DenseNDArray<std::size_t, 1>& inds) {
    if(inds.size() != dims)
      throw;
    
    std::size_t flat_ind = std::inner_product(inds.cbegin(), inds.cend(), m_strides.begin(), 0);
    return m_data[flat_ind];
  }

  bool operator==(const DenseNDArray<T, dims>& rhs) {
    return rhs.m_data == m_data;
  }
  
  auto begin() {return m_data.begin();}
  auto cbegin() {return m_data.cbegin();}
  auto end() {return m_data.end();}
  auto cend() {return m_data.cend();}
  const std::size_t size() const requires(dims == 1) {return m_data.size();}

  friend DenseNDArray<T, dims> operator+(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::plus<T>());
  }

  friend DenseNDArray<T, dims> operator+(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto plus_rhs = [&](const T& el){return std::plus<T>()(el, rhs);};
    return operator_unary(lhs, plus_rhs);
  }

  friend DenseNDArray<T, dims> operator-(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::minus<T>());
  }

  friend DenseNDArray<T, dims> operator-(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto minus_rhs = [&](const T& el){return std::minus<T>()(el, rhs);};
    return operator_unary(lhs, minus_rhs);
  }

  friend DenseNDArray<T, dims> operator*(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::multiplies<T>());
  }

  friend DenseNDArray<T, dims> operator*(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto multiplies_rhs = [&](const T& el){return std::multiplies<T>()(el, rhs);};
    return operator_unary(lhs, multiplies_rhs);
  }

  friend DenseNDArray<T, dims> operator/(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::divides<T>());
  }

  friend DenseNDArray<T, dims> operator/(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto divides_rhs = [&](const T& el){return std::divides<T>()(el, rhs);};
    return operator_unary(lhs, divides_rhs);
  }

private:

  stride_t m_strides = {};
  data_t m_data = {};

  friend DenseNDArray<T, dims> operator_binary(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs,
					       std::function<T(const T&, const T&)> binary_op) {
    DenseNDArray<T, dims> result(lhs.m_shape);
    std::transform(lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(), binary_op);
    return result;
  }

  friend DenseNDArray<T, dims> operator_unary(const DenseNDArray<T, dims>& arg, std::function<T(const T&)> unary_op) {
    DenseNDArray<T, dims> result(arg.m_shape);
    std::transform(arg.m_data.begin(), arg.m_data.end(), result.m_data.begin(), unary_op);
    return result;
  }
};

// Some type shortcuts
template <class T>
using DenseVector = DenseNDArray<T, 1>;

using IndexVector = DenseVector<std::size_t>;

#endif
