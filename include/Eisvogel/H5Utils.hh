#ifndef __H5UTILS_HH
#define __H5UTILS_HH

#include "hdf5.h"
#include <string>
#include <filesystem>

namespace H5Utils {

  hid_t open_or_create_file(std::string filepath);
  herr_t close_file(hid_t file_id);
  herr_t make_and_write_dataset(hid_t file_id, const char* name, int rank, const hsize_t* shape, hid_t dtype, const void* buffer);
  herr_t read_dataset(hid_t file_id, const char* name, hid_t dtype, void* buffer);
}

#endif
