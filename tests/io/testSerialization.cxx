#include "Serialization.hh"
#include "Eisvogel/Common.hh"
#include "DenseNDArray.hh"
#include "SparseNDArray.hh"

#include <fstream>
#include <iostream>
#include <chrono>
#include <random>

int test_serialization_vector(std::string ser_path, std::size_t length) {

  using vector_t = std::vector<scalar_t>;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(1.0, 2.0);
  
  // Generate random vector
  vector_t vec;
  for(std::size_t ind = 0; ind < length; ind++) {
    vec.push_back(dis(gen));
  }

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  
  std::fstream ofs;
  ofs.open(ser_path, std::ios::out | std::ios::binary);  
  stor::DefaultSerializer oser;
  oser.serialize(ofs, vec);
  ofs.close();

  std::fstream ifs;
  ifs.open(ser_path, std::ios::in | std::ios::binary);
  stor::DefaultSerializer iser; 
  vector_t res = iser.deserialize<vector_t>(ifs);
  ifs.close();  

  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> time_span = duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "std::vector --> Completed in " << time_span.count() << " seconds." << std::endl;
  
  for(std::size_t ind = 0; ind < length; ind++) {
    if(vec[ind] != res[ind]) {
      throw std::runtime_error("Error: mistake");
    }
  }

  return 0;
}

template <std::size_t dims>
int test_serialization_sparse_array(std::string ser_path, std::size_t size) {
  
  // Fill random sparse array
  std::array<std::size_t, dims> shape;
  for(std::size_t dim = 0; dim < dims; dim++) {
    shape[dim] = size;
  }
  DenseNDArray<scalar_t, dims> darr(shape, 1.0);

  auto to_keep = [](float value) -> bool {
    return value != 0.0;
  };
  SparseNDArray<scalar_t, dims> sparr = SparseNDArray<scalar_t, dims>::From(darr, to_keep, 0.0);

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  
  std::fstream ofs;
  ofs.open(ser_path, std::ios::out | std::ios::binary);  
  stor::DefaultSerializer oser;
  oser.serialize(ofs, sparr);
  ofs.close();
  
  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "SparseNDArray --> Serialization completed in " << time_span.count() << " seconds." << std::endl;

  t_start = std::chrono::high_resolution_clock::now();
  
  std::fstream ifs;
  ifs.open(ser_path, std::ios::in | std::ios::binary);
  stor::DefaultSerializer iser; 
  SparseNDArray<scalar_t, dims> sparr_read = iser.deserialize<SparseNDArray<scalar_t, dims>>(ifs);
  ifs.close();  

  t_end = std::chrono::high_resolution_clock::now();
  time_span = duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "SparseNDArray --> Deserialization completed in " << time_span.count() << " seconds." << std::endl;

  DenseNDArray<scalar_t, dims> darr_read = DenseNDArray<scalar_t, dims>::From(sparr_read);

  IndexVector start_inds(dims, 0);
  IndexVector end_inds = darr_read.shape();
  for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
    IndexVector cur_ind = cnt.index();
    if(darr(cur_ind) != darr_read(cur_ind)) {
      throw std::runtime_error("Error: mistake");
    }
  }

  std::cout << "Test passed" << std::endl;
  
  return 0;
}

int main(int argc, char* argv[]) {
       
  std::string ser_path = "/tmp/ser_test.bin";
  
  // test_serialization_vector(ser_path, 1e6);
  test_serialization_sparse_array<3>(ser_path, 200);
}
