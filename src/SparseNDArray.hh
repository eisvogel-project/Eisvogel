#ifndef __SPARSE_NDARRAY_HH
#define __SPARSE_NDARRAY_HH

#include <array>
#include <map>
#include "NDArray.hh"
#include "DenseNDArray.hh"
#include "Eisvogel/IteratorUtils.hh"

template <class T, std::size_t dims>
class SparseNDArray : public NDArray<T, dims> {

private:
  friend struct stor::Traits<SparseNDArray<T, dims>>;
  friend DenseNDArray<T, dims> DenseNDArray<T, dims>::From(const SparseNDArray<T, dims>& sparse_arr);
  
public:
  using shape_t = typename NDArray<T, dims>::shape_t;
  using index_t = std::array<std::size_t, dims>;
  using data_t = std::map<index_t, T>;

  SparseNDArray(const shape_t& shape, const T& default_value) : NDArray<T, dims>(shape), m_default_value(default_value) { }

  template<typename KeeperT>
  static SparseNDArray<T, dims> From(const DenseNDArray<T, dims>& dense_arr, KeeperT to_keep, const T& default_value) {

    SparseNDArray<T, dims> retval(dense_arr.shape(), default_value);

    std::size_t num_entries = 0;
    
    IndexVector start_inds(dims, 0);
    IndexVector end_inds = dense_arr.shape();
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      IndexVector cur_ind = cnt.index();
      T cur_val = dense_arr(cur_ind);

      if(to_keep(cur_val)) {
	retval.set(cur_ind, cur_val);
	num_entries++;
      }
    }

    std::cout << "kept " << num_entries << " in sparse array" << std::endl;
    
    return retval;    
  }

  std::size_t NumEntries() {
    return m_data.size();
  }
  
  T& operator()(IndexVector& inds) {

    // TODO: this will be changed when the fixed-dim vectors are available
    // Can then simply access the internal data of the fixed-dim vector class instead
    // of having to build a new index all the time like this
    index_t key;
    for(std::size_t dim = 0; dim < dims; dim++) {
      key[dim] = inds(dim);
    }
    
    auto it = m_data.find(key);
    if(it == m_data.end()) {
        return m_default_value;
    }
    return it -> second;
  }

private:

  void set(IndexVector& inds, T value) {
    
    index_t key;
    for(std::size_t dim = 0; dim < dims; dim++) {
      key[dim] = inds(dim);
    }
    
    m_data[key] = value;
  }
  
  T m_default_value;
  data_t m_data = {};  
};

namespace stor {

  template<typename T, std::size_t dims>
  struct Traits<SparseNDArray<T, dims>> {
    using type = SparseNDArray<T, dims>;
    using shape_t = typename type::shape_t;
    using data_t = typename type::data_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<T>::serialize(stream, val.m_default_value);
      Traits<shape_t>::serialize(stream, val.m_shape);
      Traits<data_t>::serialize(stream, val.m_data);
    }

    static type deserialize(std::iostream& stream) {
      T default_value = Traits<T>::deserialize(stream);
      shape_t shape = Traits<shape_t>::deserialize(stream);
      data_t data = Traits<data_t>::deserialize(stream);
      SparseNDArray<T, dims> retval(shape, default_value);
      retval.m_data = data;
      return retval;
    }
  };
}

// Some type shortcuts
template <class T>
using SparseScalarField3D = SparseNDArray<T, 3>;

#endif
