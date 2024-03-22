namespace stor {

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void DenseNDVecArraySerializer<T, dims, vec_dims>::serialize(std::fstream& stream, const type& val) {
    
    Traits<shape_t>::serialize(stream, val.m_shape);
    Traits<stride_t>::serialize(stream, val.m_strides);
    Traits<std::size_t>::serialize(stream, val.m_offset);
    Traits<data_t>::serialize(stream, *val.m_data);      
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  DenseNDVecArraySerializer<T, dims, vec_dims>::type DenseNDVecArraySerializer<T, dims, vec_dims>::deserialize(std::fstream& stream) {
    
    shape_t shape = Traits<shape_t>::deserialize(stream);
    stride_t strides = Traits<stride_t>::deserialize(stream);
    std::size_t offset = Traits<std::size_t>::deserialize(stream);
    data_t data = Traits<data_t>::deserialize(stream);
    return type(shape, strides, offset, std::move(data));
  }   
}
