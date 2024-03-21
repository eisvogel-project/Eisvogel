#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>

#include <span>

#include "Vector.hh"
#include "DenseNDVecArray.hh"

#include <stdlib.h>

int main(int argc, char* argv[]) {  

  std::shared_ptr<std::vector<float>> vec = std::make_shared<std::vector<float>>(10);

  std::span<float, 3> view(vec -> begin() + 2, 3);
  view[0] = 100;
  view[1] = 101;
  view[2] = 102;
    
  Vector<float, 3> vec1(view);

  std::cout << "vec1 = " << std::endl;
  std::cout << vec1[0] << std::endl;
  std::cout << vec1[1] << std::endl;
  std::cout << vec1[2] << std::endl;

  DenseNDVecArray<float, 2, 3> arr1({2u, 2u}, 0.0f);
  Vector<std::size_t, 2> ind{1u, 1u};  
  arr1[ind] = vec1;

  std::filesystem::path testpath = "./testVector.bin";   
  std::fstream ofs;
  ofs.open(testpath, std::ios::out | std::ios::binary);  
  stor::DefaultSerializer oser;
  oser.serialize(ofs, arr1);
  ofs.close();

  std::fstream ifs;
  ifs.open(testpath, std::ios::in | std::ios::binary);
  stor::DefaultSerializer iser;
  DenseNDVecArray<float, 2, 3> arr1_read = iser.deserialize<DenseNDVecArray<float, 2, 3>>(ifs);
  
  auto value = arr1_read[ind];
  std::cout << "retrieved value =" << std::endl;
  for(auto cur: value) {
    std::cout << cur << std::endl;
  }

  std::cout << "empty value =" << std::endl;
  for(auto cur: arr1_read[{0u, 0u}]) {
    std::cout << cur << std::endl;
  }  
  
  std::cout << "done" << std::endl;
}
