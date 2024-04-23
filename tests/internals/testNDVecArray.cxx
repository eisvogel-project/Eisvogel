#include <iostream>
#include <filesystem>
#include "NDVecArray.hh"
#include "Vector.hh"
#include "NDVecArrayStreamer.hh"
#include "IteratorUtils.hh"

int main(void) {

  Vector<float, 3> vec1(0);

  Vector<float, 3> vec2(170);

  vec1 = vec2;
  
  std::cout << vec1[0] << ", " << vec1[1] << ", " << vec1[2] << std::endl;
  
  // ----
    
  NDVecArray<float, 3, 2> arr1({2u, 2u, 2u}, 1.0f);

  arr1[{0u, 0u, 0u}] = {0.0f, 0.0f};
  arr1[{0u, 0u, 1u}] = {1.0f, -1.0f};
  arr1[{0u, 1u, 0u}] = {2.0f, -2.0f};
  arr1[{0u, 1u, 1u}] = {3.0f, -3.0f};
  arr1[{1u, 0u, 0u}] = {4.0f, -4.0f};
  arr1[{1u, 0u, 1u}] = {5.0f, -5.0f};
  arr1[{1u, 1u, 0u}] = {6.0f, -6.0f};
  arr1[{1u, 1u, 1u}] = {7.0f, -7.0f};

  std::filesystem::path testpath = "./export_from_eisvogel.bin";
  
  std::fstream iofs;
  iofs.open(testpath, std::ios::out | std::ios::binary);    
  stor::Traits<NDVecArray<float, 3, 2>>::serialize_to_numpy(iofs, arr1);
  iofs.close();
  
  // NDVecArray<float, 3, 2> arr2({3u, 3u, 3u}, 1.0f);
  
  // std::cout << "original array" << std::endl;
  // arr2.loop_over_elements(print_element);

  // arr2 = arr1;

  // std::cout << "after copy assignment" << std::endl;
  // arr2.loop_over_elements(print_element);
  
  // std::cout << "----" << std::endl;

  // // ----

  // NDVecArray<float, 3, 2> apparr_1({2u, 2u, 2u}, 1.0f);
  // NDVecArray<float, 3, 2> apparr_2({2u, 2u, 2u}, 2.0f);

  // std::cout << "before appending" << std::endl;

  // apparr_1.loop_over_elements(print_element);
  
  // apparr_1.Append<0>(apparr_2);

  // std::cout << "after appending" << std::endl;

  // apparr_1.loop_over_elements(print_element);

  // apparr_1.fill_from(arr1, {0u, 0u, 0u}, {1u, 1u, 1u}, {1u, 0u, 0u});
  
  // std::cout << "after fill_from" << std::endl;

  // apparr_1.loop_over_elements(print_element);
  
  // ----

  // for(int i = 0; i < 10; i++) {
  //   NDVecArray<float, 3, 4> largearr({400u, 400u, 400u}, 1230.0);
  //   std::cout << i << std::endl;
  // }
}
