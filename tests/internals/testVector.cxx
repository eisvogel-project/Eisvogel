#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>
#include <algorithm>

#include <span>

#include "Vector.hh"
#include "DenseNDVecArray.hh"
#include "DenseNDVecArrayStreamer.hh"
#include "DistributedNDVecArray.hh"
#include "Eisvogel/IteratorUtils.hh"

#include <stdlib.h>

int main(int argc, char* argv[]) {  

  // Vector<std::size_t, 2> start{0u, 0u};
  // Vector<std::size_t, 2> end{10u, 10u};
  // Vector<std::size_t, 2> chunk_size{100u, 300u};
  // auto print_chunk = [](Vector<std::size_t, 2>& chunk_start, Vector<std::size_t, 2>& chunk_end) {
  //   std::cout << "--------" << std::endl;
  //   std::cout << "chunk_start = " << chunk_start[0] << ", " << chunk_start[1] << std::endl;
  //   std::cout << "chunk_end = " << chunk_end[0] << ", " << chunk_end[1] << std::endl;
  //   std::cout << "--------" << std::endl;
  // };

  // auto printer = [](Vector<std::size_t, 2>& cur) {
  //   std::cout << cur[0] << ", " << cur[1]  << std::endl;
  // };
  
  // loop_over_chunks(start, end, chunk_size, print_chunk);
  // loop_over_elements(start, end, printer);
  
  Vector<float, 2> vec1{1.f, 2.f};
  
  // 3-dim array of 2-dim (field) vectors
  DenseNDVecArray<float, 3, 2> arr1({400u, 400u, 400u}, 0.0f);  
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

  Vector<std::size_t, 3> start{0u, 3u, 0u};
  Vector<std::size_t, 3> end{3u, 4u, 3u};
  auto printer = [](Vector<std::size_t, 3>& cur) {
    std::cout << cur[0] << ", " << cur[1] << ", " << cur[2] << std::endl;
  };
  loop_over_elements(start, end, printer);
  
  std::filesystem::path testpath = "./testVector.bin";   
  std::fstream ofs;
  ofs.open(testpath, std::ios::out | std::ios::binary);  
  //stor::DenseNDVecArrayStreamer<float, 3, 2>::serialize(ofs, arr1);
  stor::DenseNDVecArrayStreamer<float, 3, 2>::serialize_suppress_zero(ofs, arr1);
  ofs.close();

  std::fstream ifs;
  ifs.open(testpath, std::ios::in | std::ios::binary);
  DenseNDVecArray<float, 3, 2> arr1_read = stor::DenseNDVecArrayStreamer<float, 3, 2>::deserialize_suppress_zero(ifs);
  
  auto value = arr1_read[{350u, 350u, 1u}];
  std::cout << "retrieved value =" << std::endl;
  for(auto cur: value) {
    std::cout << cur << std::endl;
  }
    
  std::cout << "done" << std::endl;
}
