template <>
struct Traits<uint32_t> {
  using type = uint32_t;
  using ser_type = type;
  
  static void serialize(std::iostream& stream, const type& val) {
    ser_type outval = htonl(val);
    stream.write((char*)&outval, sizeof(outval));
  }
  
  static type deserialize(std::iostream& stream) {
    ser_type inval;
    stream.read((char*)&inval, sizeof(inval));
    return ntohl(inval);
  }
};

template <>
struct Traits<uint64_t> {
  using type = uint64_t;
  using ser_type = type;
  
  static void serialize(std::iostream& stream, const type& val) {
    ser_type outval = htonll(val);
    stream.write((char*)&outval, sizeof(outval));
  }
  
  static type deserialize(std::iostream& stream) {
    ser_type inval;
    stream.read((char*)&inval, sizeof(inval));
    return ntohll(inval);
  }
};

template <>
struct Traits<float> {
  using type = float;
  using ser_type = uint32_t;
  
  static void serialize(std::iostream& stream, const type& val) {
    const ser_type ser_val = reinterpret_cast<const ser_type&>(val);
    Traits<ser_type>::serialize(stream, ser_val);
  }
  
  static type deserialize(std::iostream& stream) {
    ser_type ser_val = Traits<ser_type>::deserialize(stream);
    type retval = 0.0f;
    std::memcpy(&retval, &ser_val, sizeof(ser_val));
    return retval;
  }
};

template<>
struct Traits<double> {
  using type = double;
  using ser_type = uint64_t;
  
  static void serialize(std::iostream& stream, const type& val) {
    const ser_type ser_val = reinterpret_cast<const ser_type&>(val);
    Traits<ser_type>::serialize(stream, ser_val);
  }
  
  static type deserialize(std::iostream& stream) {
    ser_type ser_val = Traits<ser_type>::deserialize(stream);
    type retval = 0.0;
    std::memcpy(&retval, &ser_val, sizeof(ser_val));
    return retval;
  }
};

template<>
struct Traits<std::vector<float>> {
  using type = std::vector<float>;
  using ser_type = uint32_t;
  
  static void serialize(std::iostream& stream, const type& val, std::size_t block_size = 10000) {
    std::size_t vec_size = val.size();
    Traits<std::size_t>::serialize(stream, vec_size);
    std::size_t vec_ind = 0;
    while(vec_ind < vec_size) {
      ser_type outbuf[std::min(vec_size - vec_ind, block_size)];
      for(std::size_t buf_ind = 0; (buf_ind < block_size) && (vec_ind < vec_size); buf_ind++, vec_ind++) {
	ser_type ser_val = reinterpret_cast<const ser_type&>(val[vec_ind]);
	outbuf[buf_ind] = htonl(ser_val);
      }
      stream.write((char*)&outbuf, sizeof(outbuf));
    }
  }
  
  static type deserialize(std::iostream& stream, std::size_t block_size = 10000) {
    std::size_t vec_size = Traits<std::size_t>::deserialize(stream);
    std::size_t vec_ind = 0;
    type val;
    while(vec_ind < vec_size) {
      ser_type inbuf[std::min(vec_size - vec_ind, block_size)];
      stream.read((char*)&inbuf, sizeof(inbuf));
      for(std::size_t buf_ind = 0; (buf_ind < block_size) && (vec_ind < vec_size); buf_ind++, vec_ind++) {
	ser_type ser_val = ntohl(inbuf[buf_ind]);
	float cur_val = 0.0f;
	std::memcpy(&cur_val, &ser_val, sizeof(ser_val));
	val.push_back(cur_val);
      }
    }
    return val;
  }
};

template<std::size_t n>
struct Traits<std::array<float, n>> {
  
};

// For general vectors
template <typename T>
struct Traits<std::vector<T>> {
  using type = std::vector<T>;
  
  static void serialize(std::iostream& stream, const type& val) {
    Traits<std::size_t>::serialize(stream, val.size());
    for(const T& cur : val) {
      Traits<T>::serialize(stream, cur);
    }
  }
  
  static type deserialize(std::iostream& stream) {
    type val;
    std::size_t size = Traits<std::size_t>::deserialize(stream);
    // TODO: speed this up by deserializing the entire array instead of values one by one
    for(std::size_t ind = 0; ind < size; ind++) {
      val.push_back(Traits<T>::deserialize(stream));
    }
    return val;
  }
};

// For general arrays
template <typename T, std::size_t n>
struct Traits<std::array<T, n>> {
  using type = std::array<T, n>;
  
  static void serialize(std::iostream& stream, const type& val) {
    // TODO: speed this up by serializing the entire array instead of values one by one
    for(const T& cur : val) {
      Traits<T>::serialize(stream, cur);
    }
  }
  
  static type deserialize(std::iostream& stream) {
    type val;
    // TODO: speed this up by serializing the entire array instead of values one by one
    for(std::size_t ind = 0; ind < n; ind++) {
      val[ind] = Traits<T>::deserialize(stream);
    }
    return val;
  }
};

template <>
struct Traits<std::string> {
  using type = std::string;
  
  static void serialize(std::iostream& stream, const type& val) {
    std::size_t num_chars = val.size();
    Traits<std::size_t>::serialize(stream, num_chars);
    stream.write(val.data(), num_chars);
  }
  
  static type deserialize(std::iostream& stream) {
    std::size_t num_chars = Traits<std::size_t>::deserialize(stream);
    std::vector<char> string_data(num_chars);
    stream.read(string_data.data(), num_chars);
    return std::string(string_data.begin(), string_data.end());
  }
};
