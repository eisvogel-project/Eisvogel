#include <iostream>
#include "NDArray.hh"
#include "DenseWeightingFieldProvider.hh"

int main(void) {

  unsigned int len_t = 10;

  int n = 10;

  NDArray<float, 2, 2> test(0.0);

  for(std::size_t i = 0; i < 2; i++) {
    for(std::size_t j = 0; j < 2; j++) {
      test(i, j) = i + j;
    }
  }

  NDArray<float, 2, 2> test2(1.0);

  //NDArray<float, 2, 2> res = test + test2;
  NDArray<float, 2, 2> res = test + 1;

  std::cout << res(0, 0) << std::endl;
  std::cout << res(0, 1) << std::endl;
  std::cout << res(1, 0) << std::endl;
  std::cout << res(1, 1) << std::endl;

  return 0;
}
