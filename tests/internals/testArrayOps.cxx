#include <iostream>
#include "NDArrayOperations.hh"

int main(int argc, char* argv[]) {

  DenseNDArray<float, 2> darr_1({3, 3}, -1.0);
  DenseNDArray<float, 2> darr_2({6, 3}, 1.0);

  std::vector<DenseNDArray<float, 2>> arrs;
  arrs.push_back(darr_1);
  arrs.push_back(darr_2);
  arrs.push_back(darr_2);
  
  DenseNDArray<float, 2> darr_res_1 = NDArrayOps::concatenate(arrs, 0);

  darr_res_1.print();
  
  std::cout << "concatenated shape" << std::endl;
  for(auto& cur : darr_res_1.shape()) {
    std::cout << cur << std::endl;
  }
  
  std::cout << "done" << std::endl;
}
