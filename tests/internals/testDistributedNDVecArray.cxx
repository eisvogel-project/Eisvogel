#include <iostream>
#include <fstream>
#include <tuple>

#include "NDVecArray.hh"
#include "DistributedNDVecArray.hh"

void fill_distributed_array(std::filesystem::path workdir) {

  Vector<std::size_t, 3> buffer_shape{400u, 400u, 400u};
  NDVecArray<float, 3, 2> chunk_buffer(buffer_shape, 1.0);

  // create a new distributed array
  Vector<std::size_t, 3> streamer_chunk_size{1u, stor::INFTY, stor::INFTY};
  DistributedNDVecArray<NDVecArray, float, 3, 2> darr(workdir, 10, buffer_shape, streamer_chunk_size);

  Vector<std::size_t, 3> start_ind{0u, 0u, 0u};
  darr.RegisterChunk(start_ind, chunk_buffer);
  
}

void read_distributed_array(std::filesystem::path workdir) {

  using darr_type = DistributedNDVecArray<NDVecArray, float, 3, 2>;
  
  darr_type darr(workdir);

  darr_type::view_t elem = darr[{1u, 2u, 3u}];

  std::cout << elem[0] << ", " << elem[1] << std::endl;
}

int main(int argc, char* argv[]) {

  auto print_element = [](VectorView<float, 2> elem) {
    std::cout << elem[0] << ", " << elem[1] << std::endl;
  };

  std::filesystem::path workdir = "./darr_test";
  if(std::filesystem::exists(workdir)) {
    std::filesystem::remove_all(workdir);
  }
  
  fill_distributed_array(workdir);
  read_distributed_array(workdir);
  
  std::cout << "done" << std::endl;
}
