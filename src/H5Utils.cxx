#include "Eisvogel/H5Utils.hh"

namespace H5Utils {

  hid_t open_or_create_file(std::string filepath) {
    hid_t file_id;
    
    if(std::filesystem::exists(filepath)) {
      file_id = H5Fopen(filepath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    } else {
      file_id = H5Fcreate(filepath.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    }
    
    return file_id;
  }

  herr_t close_file(hid_t file_id) {
    return H5Fclose(file_id);
  }
  
  herr_t make_and_write_dataset(hid_t file_id, const char* name, int rank, const hsize_t* shape, hid_t dtype, const void* buffer) {
    hid_t dataspace_id = H5Screate_simple(rank, shape, NULL);
    hid_t dataset_id = H5Dcreate2(file_id, name, dtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Dwrite(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    return status;
  }

  herr_t read_dataset(hid_t file_id, const char* name, hid_t dtype, void* buffer) {
    hid_t dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    hid_t dataspace_id = H5Dget_space(dataset_id);   
    herr_t status = H5Dread(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    return status;
  }
}
