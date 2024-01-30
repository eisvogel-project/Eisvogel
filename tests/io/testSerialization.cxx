#include "Eisvogel/Serialization.hh"
#include "Eisvogel/Common.hh"

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
  stor::Serializer oser(ofs);
  oser.serialize(vec);
  ofs.close();

  std::fstream ifs;
  ifs.open(ser_path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);    
  vector_t res = iser.deserialize<vector_t>();
  ifs.close();  

  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> time_span = duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "Completed in " << time_span.count() << " seconds." << std::endl;
  
  for(std::size_t ind = 0; ind < length; ind++) {
    if(vec[ind] != res[ind]) {
      throw std::runtime_error("Error: mistake");
    }
  }

  return 0;
}

int main(int argc, char* argv[]) {
       
  std::string ser_path = "ser_test.bin";
  test_serialization_vector(ser_path, 1e6);
}
