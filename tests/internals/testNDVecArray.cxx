#include <iostream>
#include "NDVecArray.hh"
#include "Vector.hh"

#include "Eisvogel/IteratorUtils.hh"

int main(int argc, char* argv[]) {

  Vector<float, 3> vec1(0);

  Vector<float, 3> vec2(170);

  vec1 = vec2;
  
  std::cout << vec1[0] << ", " << vec1[1] << ", " << vec1[2] << std::endl;
  return 1;
  
  // ----
  
  auto print_element = [](VectorView<float, 2> elem) {
    std::cout << elem[0] << ", " << elem[1] << std::endl;
  };  
  
  NDVecArray<float, 3, 2> arr1({2u, 2u, 2u}, 1.0f);

  arr1[{0u, 0u, 0u}] = {0.0f, 0.0f};
  arr1[{0u, 0u, 1u}] = {1.0f, -1.0f};
  arr1[{0u, 1u, 0u}] = {2.0f, -2.0f};
  arr1[{0u, 1u, 1u}] = {3.0f, -3.0f};
  arr1[{1u, 0u, 0u}] = {4.0f, -4.0f};
  arr1[{1u, 0u, 1u}] = {5.0f, -5.0f};
  arr1[{1u, 1u, 0u}] = {6.0f, -6.0f};
  arr1[{1u, 1u, 1u}] = {7.0f, -7.0f};

  NDVecArray<float, 3, 2> arr2({3u, 3u, 3u}, 1.0f);
  
  std::cout << "original array" << std::endl;
  arr2.loop_over_elements(print_element);

  arr2 = arr1;

  std::cout << "after copy assignment" << std::endl;
  arr2.loop_over_elements(print_element);
  
  std::cout << "----" << std::endl;
  
}
