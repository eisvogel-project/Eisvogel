#include <iostream>
#include <string>
#include "Eisvogel/H5Serialization.hh"
#include "Eisvogel/NDArray.hh"

int main(void) {

  std::string outpath = "test.h5";

  {
    h5stor::H5Serializer ser(outpath);
    
    // Test (de)serialization of std::array
    std::array<float, 3> testarr = {1, 2, 3};
    ser.serialize(testarr, "testarr");
  }

  {
    h5stor::H5Serializer ser(outpath);
    std::array<float, 3> read_arr = ser.deserialize<std::array<float, 3>>("testarr");
    for(auto cur: read_arr) {
      std::cout << cur << std::endl;
    }
  }
    
  // // Build a simple array for testing purposes
  // DenseNDArray<float, 2> testarr({2, 2}, 1.0);

  // std::cout << testarr(1,1) << std::endl;

  // ser.serialize(testarr, "testarr");

  // ser.deserialize<DenseNDArray<float, 2>>("testarr");
  
  std::cout << "Done" << std::endl;
  return 0;
}
