#include <iostream>
#include <string>
#include "Eisvogel/H5Serialization.hh"
#include "Eisvogel/NDArray.hh"

int main(void) {

  std::string outpath = "test.h5";
  h5stor::H5Serializer ser(outpath);

  // Build a simple array for testing purposes
  DenseNDArray<float, 2> testarr({2, 2}, 1.0);

  std::cout << testarr(1,1) << std::endl;

  ser.serialize(testarr, "testarr");
  
  std::cout << "Done" << std::endl;
  return 0;
}
