#include <iostream>
#include "NDVecArray.hh"

#include "Eisvogel/IteratorUtils.hh"

int main(int argc, char* argv[]) {

  NDVecArray<float, 3, 2> arr1({2u, 2u, 2u}, 1.0f);

  arr1[{0u, 0u, 0u}] = {0.0f, 0.0f};
  arr1[{0u, 0u, 1u}] = {1.0f, -1.0f};
  arr1[{0u, 1u, 0u}] = {2.0f, -2.0f};
  arr1[{0u, 1u, 1u}] = {3.0f, -3.0f};
  arr1[{1u, 0u, 0u}] = {4.0f, -4.0f};
  arr1[{1u, 0u, 1u}] = {5.0f, -5.0f};
  arr1[{1u, 1u, 0u}] = {6.0f, -6.0f};
  arr1[{1u, 1u, 1u}] = {7.0f, -7.0f};

  auto print_element = [](VectorView<float, 2> elem) {
    std::cout << elem[0] << ", " << elem[1] << std::endl;
  };  
  arr1.loop_over_elements(print_element);

  std::cout << "----" << std::endl;
  
  auto print_element_by_index = [&](Vector<std::size_t, 3> ind) {
    std::cout << arr1[ind][0] << ", " << arr1[ind][1] << std::endl;
  };
  index_loop_over_array_elements(arr1, print_element_by_index);
}
