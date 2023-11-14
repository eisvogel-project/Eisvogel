#include <iostream>
#include <string>
#include "Eisvogel/H5Serialization.hh"

int main(void) {

  std::string outpath = "test.h5";
  h5stor::H5Serializer ser(outpath);
  
  std::cout << "Done" << std::endl;
  return 0;
}
