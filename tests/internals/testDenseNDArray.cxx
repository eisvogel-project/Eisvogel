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

  DenseNDArray<float, 2> darr_restricted = NDArrayOps::range(darr_res, {3, 0}, {9, 3});
  std::cout << "restricted" << std::endl;

  for(auto cur : darr_restricted.shape()) {
    std::cout << cur << std::endl;
  }

  std::cout << "vals" << std::endl;
  std::cout << darr_restricted(5, 1) << std::endl;
  
  std::cout << "done" << std::endl;
}
