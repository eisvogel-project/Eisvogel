#ifndef __DENSE_NDARRAY_HH
#define __DENSE_NDARRAY_HH

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include "NDArray.hh"
#include "Serialization.hh"

template <class T, std::size_t dims>
class SparseNDArray;

// TODO: probably also want a fixed-sized version of NDArray that allows to specify the individual dimensions
// Useful for e.g. 3-dim vectors that don't need to be dynamic

template <class T, std::size_t dims>
class DenseNDArray : public NDArray<T, dims> {

public:
  using shape_t = typename NDArray<T, dims>::shape_t;

private:
  using stride_t = std::array<std::size_t, dims + 1>;
  using data_t = std::vector<T>;

private:
  friend struct stor::Traits<DenseNDArray<T, dims>>;
  
  template <typename T0, std::size_t dims0, typename T1, typename T2>
  friend inline DenseNDArray<T0, dims0> operator_binary(const DenseNDArray<T1, dims0>& lhs, const DenseNDArray<T2, dims0>& rhs, auto binary_op);
  
  template <typename T0, std::size_t dims0, typename T1>
  friend inline DenseNDArray<T0, dims0> operator_unary(const DenseNDArray<T1, dims0>& arg, auto unary_op);

  template <typename T1> friend class DenseNDArray<T, dims>;

public:

  using type = T;

  DenseNDArray(std::size_t size, const T& value) requires(dims == 1) : DenseNDArray(shape_t({size}), value) { }

  // TODO: This will move to use a fixed-size vector to specify the shape (and dims)
  DenseNDArray(const shape_t& shape, const T& value) : NDArray<T, dims>(shape) {
    m_strides[0] = 1;
    std::partial_sum(shape.begin(), shape.end(), m_strides.begin() + 1, std::multiplies<std::size_t>());
    m_data.resize(m_strides.back(), value);
  }
  
  static DenseNDArray<T, dims> From(const SparseNDArray<T, dims>& sparse_arr) {
    DenseNDArray<T, dims> retval(sparse_arr.shape(), sparse_arr.m_default_value);

    for (auto const& [key, val] : sparse_arr.m_data) {
      retval(key) = val;
    }
    
    return retval;
  }

  static DenseNDArray<T, dims> FromSparseFile(std::iostream& stream) {
    T default_value = stor::Traits<T>::deserialize(stream);
    shape_t shape = stor::Traits<shape_t>::deserialize(stream);
    DenseNDArray<T, dims> retval(shape, default_value);

    std::vector<std::array<std::size_t, dims>> keys = stor::Traits<std::vector<std::array<std::size_t, dims>>>::deserialize(stream);
    std::vector<T> values = stor::Traits<std::vector<T>>::deserialize(stream);

    for(std::size_t el_ind = 0; el_ind < keys.size(); el_ind++) {
      retval(keys[el_ind]) = values[el_ind];
    }
    
    return retval;
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

  template <typename T1>
  DenseNDArray(std::vector<T1>&& data) requires(dims == 1) : 
    NDArray<T, 1>({data.size()}), m_strides({1, data.size()}), m_data(data.begin(), data.end()) { }

  template <typename T1, std::size_t length>
  DenseNDArray(const std::array<T1, length>& data) requires(dims == 1) : 
    NDArray<T, 1>({length}), m_strides({1, length}), m_data(data.cbegin(), data.cend()) { }

  template <typename T1>
  DenseNDArray(const DenseNDArray<T1, dims>& other) : 
    NDArray<T, dims>(other.m_shape), m_strides(other.m_strides), m_data(other.m_data.cbegin(), other.m_data.cend()) { }

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
  // TODO: move to using vectors with compile-time size here and make sure the size matches
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
  
  T& operator()(const std::array<std::size_t, dims>& inds) {
    std::size_t flat_ind = std::inner_product(inds.cbegin(), inds.cend(), m_strides.begin(), 0);
    return m_data[flat_ind];    
  }

  // indexing for special (low-dimensional) cases
  inline const T& operator()(std::size_t ind) const requires(dims == 1) {
    return m_data[ind];
  }

  inline T& operator()(std::size_t ind) requires(dims == 1) {
    return m_data[ind];
  }

  inline const T& operator()(DenseNDArray<std::size_t, 1>& inds) const requires(dims == 3) {
    return m_data[inds.m_data[0] * m_strides[0] + inds.m_data[1] * m_strides[1] + inds.m_data[2] * m_strides[2]];
  }

  inline T& operator()(DenseNDArray<std::size_t, 1>& inds) requires(dims == 3) {
    return m_data[inds.m_data[0] * m_strides[0] + inds.m_data[1] * m_strides[1] + inds.m_data[2] * m_strides[2]];
  }

  bool operator==(const DenseNDArray<T, dims>& rhs) {
    return rhs.m_data == m_data;
  }

  // to assign blocks of data in an efficient way
  void copy_from(const DenseNDArray<T, dims>& src,
		 const DenseNDArray<std::size_t, 1>& ind_src_start, const DenseNDArray<std::size_t, 1>& ind_src_stop,
		 const DenseNDArray<std::size_t, 1>& ind_dest_start, const DenseNDArray<std::size_t, 1>& ind_dest_stop) {

    
  }
  
  // template <std::size_t len>
  // bool operator==(const std::array<T, len>& rhs) const requires(dims == 1) {
  //   if(len != m_data.size()) {
  //     return false;
  //   }
  //   for(std::size_t ind = 0; ind < len; ind++) {
  //     if(m_data[ind] != rhs[ind]) {
  // 	return false;
  //     }
  //   }
  //   return true;
  // }

  // printing
  void print() const {
    for(const T& cur: m_data) {
      std::cout << cur << " ";
    }
    std::cout << std::endl;
  }
  
  auto begin() {return m_data.begin();}
  auto cbegin() {return m_data.cbegin();}
  auto begin() const {return m_data.cbegin();}
  auto end() {return m_data.end();}
  auto cend() {return m_data.cend();}
  auto end() const {return m_data.cend();}
  std::size_t size() const requires(dims == 1) {return m_data.size();}
  std::size_t volume() const {return m_strides.back();}

  friend DenseNDArray<T, dims> operator+(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary<T, dims, T, T>(lhs, rhs, std::plus<>());
  }

  friend DenseNDArray<T, dims> operator+(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto plus_rhs = [&](const T& el){return std::plus<>()(el, rhs);};
    return operator_unary<T, dims, T>(lhs, plus_rhs);
  }

  friend DenseNDArray<T, dims> operator-(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary<T, dims, T, T>(lhs, rhs, std::minus<>());
  }

  friend DenseNDArray<T, dims> operator-(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto minus_rhs = [&](const T& el){return std::minus<>()(el, rhs);};
    return operator_unary<T, dims, T>(lhs, minus_rhs);
  }

  friend DenseNDArray<T, dims> operator*(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary<T, dims, T, T>(lhs, rhs, std::multiplies<>());
  }

  friend DenseNDArray<T, dims> operator*(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto multiplies_rhs = [&](const T& el){return std::multiplies<>()(el, rhs);};
    return operator_unary<T, dims, T>(lhs, multiplies_rhs);
  }

  friend DenseNDArray<T, dims> operator/(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    return operator_binary<T, dims, T, T>(lhs, rhs, std::divides<>());
  }

  friend DenseNDArray<T, dims> operator/(const DenseNDArray<T, dims>& lhs, const T& rhs) {
    auto divides_rhs = [&](const T& el){return std::divides<>()(el, rhs);};
    return operator_unary<T, dims, T>(lhs, divides_rhs);
  }

  friend DenseNDArray<bool, dims> operator<(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    auto elementwise_lt = [&](const T& el_lhs, const T& el_rhs){return el_lhs < el_rhs;};
    return operator_binary<bool, dims, T, T>(lhs, rhs, elementwise_lt);
  }

  friend DenseNDArray<bool, dims> operator<=(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    auto elementwise_leq = [&](const T& el_lhs, const T& el_rhs){return el_lhs <= el_rhs;};
    return operator_binary<bool, dims, T, T>(lhs, rhs, elementwise_leq);
  }

  friend DenseNDArray<bool, dims> operator>(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    auto elementwise_gt = [&](const T& el_lhs, const T& el_rhs){return el_lhs > el_rhs;};
    return operator_binary<bool, dims, T, T>(lhs, rhs, elementwise_gt);
  }

  friend DenseNDArray<bool, dims> operator>=(const DenseNDArray<T, dims>& lhs, const DenseNDArray<T, dims>& rhs) {
    auto elementwise_geq = [&](const T& el_lhs, const T& el_rhs){return el_lhs >= el_rhs;};
    return operator_binary<bool, dims, T, T>(lhs, rhs, elementwise_geq);
  }
  
private:

  stride_t m_strides = {};
  data_t m_data = {};
};

template <typename T0, std::size_t dims0, typename T1, typename T2>
inline DenseNDArray<T0, dims0> operator_binary(const DenseNDArray<T1, dims0>& lhs, const DenseNDArray<T2, dims0>& rhs,
					       auto binary_op) {
  DenseNDArray<T0, dims0> result(lhs.m_shape, T0());
  std::transform(lhs.m_data.begin(), lhs.m_data.end(), rhs.m_data.begin(), result.m_data.begin(), binary_op);
  return result;
}

template <typename T0, std::size_t dims0, typename T1>
inline DenseNDArray<T0, dims0> operator_unary(const DenseNDArray<T1, dims0>& arg, auto unary_op) {
  DenseNDArray<T0, dims0> result(arg.m_shape, T0());
  std::transform(arg.m_data.begin(), arg.m_data.end(), result.m_data.begin(), unary_op);
  return result;
}

template <std::size_t dims>
bool all(const DenseNDArray<bool, dims>& mask) {
  for(bool cur : mask) {
    if(cur == false) {
      return false;
    }
  }
  return true;
}

namespace stor {

  template<typename T, std::size_t dims>
  struct Traits<DenseNDArray<T, dims>> {
    using type = DenseNDArray<T, dims>;
    using shape_t = typename type::shape_t;
    using data_t = typename type::data_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<shape_t>::serialize(stream, val.m_shape);
      Traits<data_t>::serialize(stream, val.m_data);
    }

    static type deserialize(std::iostream& stream) {
      shape_t shape = Traits<shape_t>::deserialize(stream);
      data_t data = Traits<data_t>::deserialize(stream);
      return DenseNDArray<T, dims>(shape, std::move(data));
    }
  };
}

// Some type shortcuts
template <class T>
using DenseVector = DenseNDArray<T, 1>;

using IndexVector = DenseVector<std::size_t>;
using GridVector = DenseVector<unsigned int>;

template <class T>
using ScalarField3D = DenseNDArray<T, 3>;

template <class T>
using ScalarField2D = DenseNDArray<T, 2>;

#endif
