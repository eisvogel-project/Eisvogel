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

  template <class Derived>
  class Serializer {

  public:
    
    template <typename T>
    void serialize(std::fstream& stream, const T& value) {
      static_cast<Derived*>(this) -> template serialize_impl<T>(stream, value);
    }

    template <typename T>
    T deserialize(std::fstream& stream) {
      return static_cast<Derived*>(this) -> template deserialize_impl<T>(stream);
    }    
  };
  
  class DefaultSerializer : public Serializer<DefaultSerializer> {

  public:
    DefaultSerializer() { }
    
    template <typename T>
    void serialize_impl(std::fstream& stream, const T& value);

    template <typename T>
    T deserialize_impl(std::fstream& stream);
    
  };

}

#include "Serialization.hxx"
#include "DefaultSerializationTraits.hxx"

#endif
