#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>

#include <span>

#include "Vector.hh"
#include "DenseNDVecArray.hh"
#include "DistributedNDVecArray.hh"

#include <stdlib.h>

int main(int argc, char* argv[]) {  

  Vector<float, 2> vec1{1.f, 2.f};
  
  // 3-dim array of 2-dim (field) vectors
  DenseNDVecArray<float, 3, 2> arr1({400u, 400u, 400u}, 0.0f);
  
  Vector<std::size_t, 2> ind{1u, 1u};
  arr1[{350u, 350u, 1u}] = vec1;

  std::cout << "empty value =" << std::endl;
  for(auto cur: arr1[{0u, 0u, 0u}]) {
    std::cout << cur << std::endl;
  }  

  std::cout << "filled value =" << std::endl;
  for(auto cur: arr1[{350u, 350u, 1u}]) {
    std::cout << cur << std::endl;
  }

  DenseNDVecArray<float, 3, 2> arr_view = arr1.View({350u, 350u, 1u}, {351u, 351u, 10u});

  std::cout << "empty value in view =" << std::endl;
  for(auto cur: arr_view[{0u, 0u, 3u}]) {
    std::cout << cur << std::endl;
  }
  
  std::cout << "filled value in view =" << std::endl;
  for(auto cur: arr_view[{0u, 0u, 0u}]) {
    std::cout << cur << std::endl;
  }
  
  // std::filesystem::path testpath = "./testVector.bin";   
  // std::fstream ofs;
  // ofs.open(testpath, std::ios::out | std::ios::binary);  
  // stor::DefaultSerializer oser;
  // oser.serialize(ofs, arr1);
  // ofs.close();

  // std::fstream ifs;
  // ifs.open(testpath, std::ios::in | std::ios::binary);
  // stor::DefaultSerializer iser;
  // DenseNDVecArray<float, 2, 3> arr1_read = iser.deserialize<DenseNDVecArray<float, 2, 3>>(ifs);
  
  // auto value = arr1_read[ind];
  // std::cout << "retrieved value =" << std::endl;
  // for(auto cur: value) {
  //   std::cout << cur << std::endl;
  // }
  
  std::cout << "done" << std::endl;
}
