#include <iostream>
#include "DenseNDArray.hh"
#include "NDArrayOperations.hh"

int test_concatenation() {
  
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

  return 0;
}

int test_slicing() {

  DenseNDArray<float, 2> darr_1({30, 30}, -1.0);
  DenseNDArray<float, 2> darr_2({3, 3}, 1.0);

  darr_1.copy_from(darr_2,
		   {0, 0}, {2, 2},   // index range in src array
		   {0, 0}, {2, 2});  // index range in destination array

  std::cout << darr_1(0, 0) << std::endl;
  std::cout << darr_1(0, 1) << std::endl;
  std::cout << darr_1(1, 1) << std::endl;
  
  return 0;
}

int main(int argc, char* argv[]) {

  // test_concatenation();

  test_slicing();
  
  std::cout << "done" << std::endl;
}
