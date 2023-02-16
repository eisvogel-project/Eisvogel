#ifndef __NDARRAY_HH
#define __NDARRAY_HH

#include <vector>
#include <array>
#include <numeric>
#include <algorithm>

template <class T, std::size_t dims>
class NDArray {

private:

  using shape_t = std::array<std::size_t, dims>;
  using stride_t = std::array<std::size_t, dims + 1>;
  using data_t = std::vector<T>;

public:
  
  NDArray(const shape_t& shape, const T& value = T()) : m_shape(shape) {   
    m_strides[0] = 1;
    std::partial_sum(m_shape.begin(), m_shape.end(), m_strides.begin() + 1, std::multiplies<std::size_t>());
    m_data.resize(m_strides.back(), value);
  }

  NDArray(const NDArray<T, dims>& other) : m_shape(other.m_shape), m_strides(other.m_strides), m_data(other.m_data) {}

  template <typename... Inds>
  T& operator()(Inds... inds) {
    std::size_t flat_ind = 0, dim = 0;
    (..., (flat_ind += inds * m_strides[dim++]));
    return m_data[flat_ind];
  }

  friend NDArray<T, dims> operator+(const NDArray<T, dims>& lhs, const NDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::plus<T>());
  }

  friend NDArray<T, dims> operator+(const NDArray<T, dims>& lhs, const T& rhs) {
    auto plus_rhs = [&](const T& el){return std::plus<T>()(el, rhs);};
    return operator_unary(lhs, plus_rhs);
  }

  friend NDArray<T, dims> operator-(const NDArray<T, dims>& lhs, const NDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::minus<T>());
  }

  friend NDArray<T, dims> operator-(const NDArray<T, dims>& lhs, const T& rhs) {
    auto minus_rhs = [&](const T& el){return std::minus<T>()(el, rhs);};
    return operator_unary(lhs, minus_rhs);
  }

  friend NDArray<T, dims> operator*(const NDArray<T, dims>& lhs, const NDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::multiplies<T>());
  }

  friend NDArray<T, dims> operator*(const NDArray<T, dims>& lhs, const T& rhs) {
    auto multiplies_rhs = [&](const T& el){return std::multiplies<T>()(el, rhs);};
    return operator_unary(lhs, multiplies_rhs);
  }

  friend NDArray<T, dims> operator/(const NDArray<T, dims>& lhs, const NDArray<T, dims>& rhs) {
    return operator_binary(lhs, rhs, std::divides<T>());
  }

  friend NDArray<T, dims> operator/(const NDArray<T, dims>& lhs, const T& rhs) {
    auto divides_rhs = [&](const T& el){return std::divides<T>()(el, rhs);};
    return operator_unary(lhs, divides_rhs);
  }

private:

  shape_t m_shape = {};
  stride_t m_strides = {};
  data_t m_data = {};

  friend NDArray<T, dims> operator_binary(const NDArray<T, dims>& lhs, const NDArray<T, dims>& rhs,
					  std::function<T(const T&, const T&)> binary_op) {
    NDArray<T, dims> result(lhs.m_shape);
    std::transform(lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(), binary_op);
    return result;
  }

  friend NDArray<T, dims> operator_unary(const NDArray<T, dims>& arg, std::function<T(const T&)> unary_op) {
    NDArray<T, dims> result(arg.m_shape);
    std::transform(arg.m_data.begin(), arg.m_data.end(), result.m_data.begin(), unary_op);
    return result;
  }
};

#endif
