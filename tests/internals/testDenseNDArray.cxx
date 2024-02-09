#include <iostream>
#include <fstream>
#include "Serialization.hh"
#include "DenseNDArray.hh"
#include "NDArrayOperations.hh"

int main(int argc, char* argv[]) {

  DenseNDArray<float, 2> darr_1({3, 3}, -1.0);
  DenseNDArray<float, 2> darr_2({6, 3}, 1.0);

  DenseNDArray<float, 2> darr_res = NDArrayOps::concatenate(darr_1, darr_2, 0);

  for(auto cur : darr_res.shape()) {
    std::cout << cur << std::endl;
  }

  std::cout << "vals" << std::endl;
  std::cout << darr_res(2, 2) << std::endl;
  
  std::cout << "done" << std::endl;
}
