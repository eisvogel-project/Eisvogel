#include "Eisvogel/Serialization.hh"
#include "Eisvogel/Common.hh"

#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {

  using vector_t = std::vector<scalar_t>;
  
  vector_t vec = {-1, -2, 3, 4};

  std::string ser_path = "ser_test.bin";
  
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

  for(auto cur: res) {
    std::cout << cur << std::endl;
  }
}
