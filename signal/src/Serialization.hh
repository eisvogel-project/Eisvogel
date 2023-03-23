#ifndef __SERIALIZATION_HH
#define __SERIALIZATION_HH

#include <iostream>
#include <cstdint>
#include <vector>
#include <array>
#include <arpa/inet.h>

// inspired by https://github.com/panta/seriously

namespace stor {

  // https://stackoverflow.com/a/28592202
#define htonll(x) ((1==htonl(1)) ? (x) : ((uint64_t)htonl((x) & 0xFFFFFFFF) << 32) | htonl((x) >> 32))
#define ntohll(x) ((1==ntohl(1)) ? (x) : ((uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))

  template <typename T>
  struct Traits;

  template <>
  struct Traits<uint32_t> {
    using type = uint32_t;
    using ser_type = type;

    static void serialize(std::iostream& stream, const type& val) {
      ser_type outval = htonl(val);
      stream.write((char*)&outval, sizeof(outval));
    }

    static void deserialize(std::iostream& stream, type& val) {
      ser_type inval;
      stream.read((char*)&inval, sizeof(inval));
      val = ntohl(inval);
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

    static void deserialize(std::iostream& stream, type& val) {
      ser_type inval;
      stream.read((char*)&inval, sizeof(inval));
      val = ntohll(inval);
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

    static void deserialize(std::iostream& stream, type& val) {
      ser_type ser_val;
      Traits<ser_type>::deserialize(stream, ser_val);
      val = reinterpret_cast<type&>(ser_val);
    }
  };

  template <typename T>
  struct Traits<std::vector<T>> {
    using type = std::vector<T>;
    
    static void serialize(std::iostream& stream, const type& val) {
      Traits<std::size_t>::serialize(stream, val.size());
      for(const T& cur : val) {
	Traits<T>::serialize(stream, cur);
      }
    }
    
    static void deserialize(std::iostream& stream, type& val) {
      std::size_t size;
      Traits<std::size_t>::deserialize(stream, size);
      for(std::size_t ind = 0; ind < size; ind++) {
	T cur;
	Traits<T>::deserialize(stream, cur);
	val.push_back(cur);
      }
    }
  };

  template <typename T, std::size_t n>
  struct Traits<std::array<T, n>> {
    using type = std::array<T, n>;

    static void serialize(std::iostream& stream, const type& val) {
      for(const T& cur : val) {
	Traits<T>::serialize(stream, cur);
      }
    }

    static void deserialize(std::iostream& stream, type& val) {
      for(std::size_t ind = 0; ind < n; ind++) {
	T cur;
	Traits<T>::deserialize(stream, cur);
	val[ind] = cur;
      }
    }
  };

  class Serializer {

  public:
    Serializer(std::iostream& stream) : m_stream(stream) { }
    
    template <typename T>
    friend Serializer operator<<(Serializer& ser, const T& value) {
      Traits<T>::serialize(ser.stream(), value);
      return ser;
    }

    template <typename T>
    friend Serializer operator>>(Serializer& ser, T& value) {
      Traits<T>::deserialize(ser.stream(), value);
      return ser;
    }

    std::iostream& stream() {
      return m_stream;
    }

  private:
    std::iostream& m_stream;

  };
}

#endif
