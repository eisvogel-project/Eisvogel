#ifndef __H5SERIALIZATION_HH
#define __H5SERIALIZATION_HH

#include "hdf5.h"
#include <string>

namespace h5stor {

  template <typename T>
  struct Traits;
  
  class H5Serializer {

  public:
    H5Serializer(std::string filepath) {
      m_file_id = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
