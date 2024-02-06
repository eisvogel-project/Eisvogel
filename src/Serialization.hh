#ifndef __SERIALIZATION_HH
#define __SERIALIZATION_HH

#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <array>
#include <arpa/inet.h>
#include <cstring>
#include "Eisvogel/Common.hh"

#include <iostream>

// inspired by https://github.com/panta/seriously

namespace stor {

  // https://stackoverflow.com/a/28592202
#define htonll(x) ((1==htonl(1)) ? (x) : ((uint64_t)htonl((x) & 0xFFFFFFFF) << 32) | htonl((x) >> 32))
#define ntohll(x) ((1==ntohl(1)) ? (x) : ((uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))

  template <typename T>
  struct Traits;
  
  class Serializer {

  public:
    Serializer() { }
    
    template <typename T>
    void serialize(std::fstream& stream, const T& value) {
      Traits<T>::serialize(stream, value);
    }

    template <typename T>
    T deserialize(std::fstream& stream) {
      return Traits<T>::deserialize(stream);
    }
  };

  class SparseSerializer : public Serializer {

  public:
    SparseSerializer() { }

    template <typename T>
    void serialize(std::fstream& stream, const T& value) {
      std::cout << "serializing in a sparse way!" << std::endl;
      Traits<T>::serialize(stream, value);
    }

    template <typename T>
    T deserialize(std::fstream& stream) {
      std::cout << "deserializing in a sparse way!" << std::endl;
      return Traits<T>::deserialize(stream);
    }
  };
}

#include "Serialization.hxx"

#endif
