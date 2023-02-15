#ifndef __NDARRAY_HH
#define __NDARRAY_HH

#include <valarray>
#include <functional>

template <class T, std::size_t... Shape>
class NDArray {

private:
  using stride_t = std::valarray<std::size_t>;
  using data_t = std::valarray<T>;

public:

  NDArray(const T value) {
    std::size_t seed = 1;
    strides = {seed, seed *= Shape...};

    data.resize(strides[strides.size() - 1], value);
  }

  NDArray(const NDArray<T, Shape...>& other) : strides(other.strides), data(other.data) {}

  template <typename... Inds>
  T& operator()(Inds... inds) {
    std::size_t flat_ind = 0, dim = 0;
    (..., (flat_ind += inds * strides[dim++]));

    return data[flat_ind];
  }

  friend NDArray<T, Shape...> operator+(const NDArray<T, Shape...>& lhs, const NDArray<T, Shape...>& rhs) {
    return operator_binary(lhs, rhs, std::plus<data_t>());
  }

  friend NDArray<T, Shape...> operator+(const NDArray<T, Shape...>& lhs, const T& rhs) {
    return operator_binary(lhs, rhs, std::plus<data_t>());
  }

private:

  stride_t strides;
  data_t data;

  friend NDArray<T, Shape...> operator_binary(const NDArray<T, Shape...>& lhs, const NDArray<T, Shape...>& rhs, 
					      std::function<data_t(data_t, data_t)> op) {
    NDArray<T, Shape...> result(lhs);
    result.data = op(result.data, rhs.data);
    return result;
  }
};

#endif
