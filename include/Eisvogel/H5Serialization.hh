#ifndef __H5SERIALIZATION_HH
#define __H5SERIALIZATION_HH

#include <string>
#include "hdf5.h"
#include "Eisvogel/H5Utils.hh"

namespace h5stor {

  template <typename T>
  struct Traits;

  template <std::size_t n>
  struct Traits<std::array<float, n>> {
    using type = std::array<float, n>;

    static void serialize(hid_t m_file_id, const type& val, std::string name) {
      hsize_t shape[1] = {n};
      H5Utils::make_and_write_dataset(m_file_id, name.c_str(), 1, shape, H5T_NATIVE_FLOAT, val.data());
    }

    static type deserialize(hid_t m_file_id, std::string name) {
      std::array<float, n> val;
      H5Utils::read_dataset(m_file_id, name.c_str(), H5T_NATIVE_FLOAT, val.data());
      return val;
    }
  }; 
  
  class H5Serializer {

  public:
    H5Serializer(std::string filepath) {
      m_file_id = H5Utils::open_or_create_file(filepath);
    }

    ~H5Serializer() {
      H5Utils::close_file(m_file_id);
    }
    
    template <typename T>
    void serialize(const T& value, std::string name) {
      Traits<T>::serialize(m_file_id, value, name);
    }

    template <typename T>
    T deserialize(std::string name) {
      return Traits<T>::deserialize(m_file_id, name);
    }    
    
  private:
    hid_t m_file_id;
  };
}

#endif
