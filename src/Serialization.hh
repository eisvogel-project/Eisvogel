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

// inspired by https://github.com/panta/seriously

namespace stor {

  // https://stackoverflow.com/a/28592202
#define htonll(x) ((1==htonl(1)) ? (x) : ((uint64_t)htonl((x) & 0xFFFFFFFF) << 32) | htonl((x) >> 32))
#define ntohll(x) ((1==ntohl(1)) ? (x) : ((uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))

  template <typename T>
  struct Traits;
  
  class Serializer {

  public:
    Serializer(std::fstream& stream) : m_stream(stream) { }
    
    template <typename T>
    void serialize(const T& value) {
      Traits<T>::serialize(m_stream, value);
    }

    template <typename T>
    T deserialize() {
      return Traits<T>::deserialize(m_stream);
    }

  private:
    std::fstream& m_stream;
  };
}

#include "Serialization.hxx"

#endif
