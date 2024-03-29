#include <iostream>
#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <limits>
#include <chrono>

#include <span>

#include "Vector.hh"
#include "Eisvogel/MathUtils.hh"
#include "Eisvogel/IteratorUtils.hh"
#include "NDVecArray.hh"
#include "NDVecArrayStreamer.hh"
#include "DistributedNDVecArray.hh"
#include "Eisvogel/IteratorUtils.hh"

#include <stdlib.h>

int main(int argc, char* argv[]) { 

  static constexpr std::size_t len = 400;
  
  Vector<std::size_t, 3> shape{401u, 400u, 400u};
  
  NDVecArray<float, 3, 2> arr1(shape, 150.0f);  
  std::cout << "volume = " << arr1.GetVolume() << std::endl;
  std::cout << "elements = " << arr1.GetNumberElements() << std::endl;

  Vector<std::size_t, 3> slice_shape{1u, 400u, 400u};
  
  NDVecArray<float, 3, 2> slice1 = arr1.View({400u, 0u, 0u}, {401u, 400u, 400u});

  std::cout << slice1[{0u, 0u, 0u}][0] << std::endl;
  
  Vector<float, 2> testvec{1.f, 2.f};
  arr1[{250u, 250u, 250u}] = testvec;

  std::cout << arr1[{250u, 250u, 250u}][0] << std::endl;
  std::cout << arr1[{250u, 250u, 250u}][1] << std::endl;

  NDVecArray<float, 3, 2> arr1_view = arr1.View({0u, 0u, 0u}, {400u, 400u, 400u});
  
  using in_type = float;
  using ser_type = uint32_t;

  Vector<std::size_t, 3> streamer_chunk_size{1u, stor::INFTY, stor::INFTY};
  stor::NDVecArrayStreamer<NDVecArray, float, 3, 2> streamer;
  
  std::filesystem::path testpath = "./testVector.bin";

  auto start = std::chrono::high_resolution_clock::now();
  
  {
    std::fstream iofs;
    iofs.open(testpath, std::ios::out | std::ios::binary);    
    //streamer.serialize(iofs, arr1, streamer_chunk_size, stor::StreamerMode::dense);
    streamer.serialize(iofs, arr1_view, streamer_chunk_size, stor::StreamerMode::dense);
    //streamer.serialize(iofs, slice1, streamer_chunk_size, stor::StreamerMode::dense);
    std::cout << "done writing" << std::endl;
    iofs.close();
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  // std::cout << "wrote 512MB in " << duration << std::endl;
  
  // {
  //   std::fstream iofs;
  //   iofs.open(testpath, std::ios::in | std::ios::out | std::ios::binary);    
  //   streamer.mark_as_final(iofs);
  //   std::cout << "done marking as final" << std::endl;
  //   iofs.close();
  // }
  
  {
    std::fstream iofs;
    iofs.open(testpath, std::ios::in | std::ios::out | std::ios::binary);
    streamer.append_slice<0>(iofs, slice1, stor::StreamerMode::dense);
    std::cout << "done appending slice" << std::endl;
    iofs.close();
  }

  // //Vector<std::size_t, 3> bigger_shape{4u, len, len};
  NDVecArray<float, 3, 2> arr1_read(shape, 1.0f);

  {
    std::fstream iofs;
    iofs.open(testpath, std::ios::in | std::ios::binary);    
    streamer.deserialize(iofs, arr1_read);
    std::cout << "returned" << std::endl;
    iofs.close();
  }  

  std::cout << "read array with shape = " << arr1_read.GetShape()[0] << ", " << arr1_read.GetShape()[1] << ", " << arr1_read.GetShape()[2] << std::endl;
  
  auto checker = [&](const Vector<std::size_t, 3>& ind) {
    for(std::size_t vec_ind = 0; vec_ind < 2; vec_ind++) {

      // std::cout << "Have " << arr1[ind][vec_ind] << " vs " << arr1_read[ind][vec_ind] << std::endl;
      
      if(arr1[ind][vec_ind] != arr1_read[ind][vec_ind]) {
	std::cout << "Problem at ind = " << ind[0] << ", " << ind[1] << ", " << ind[2] << " and vec_ind = " << vec_ind <<" !" << std::endl;
	std::cout << "Have " << arr1[ind][vec_ind] << " vs " << arr1_read[ind][vec_ind] << std::endl;
      }
    }   
  };

  std::cout << "checking" << std::endl;
  index_loop_over_array_elements(arr1_view, checker);

  // =====
  
  // std::size_t buflen = nullsup::calculate_required_buflen(arr1);
  // std::cout << "buflen = " << buflen << std::endl;
  // std::cout << "buflen = " << nullsup::calculate_required_buflen(shape, 2) << std::endl;
  // std::vector<ser_type> buffer(buflen, 0);

  // std::size_t elems_written = nullsup::suppress_zero(arr1, std::span<ser_type>(buffer));

  // std::cout << "wrote " << elems_written << " elements into output buffer of length " << buflen << std::endl;
  // for(std::size_t ind = 0; ind < elems_written; ind++) {
  //   std::cout << buffer[ind] << std::endl;
  // }

  // NDVecArray<float, 3, 2> arr1_read({400u, 400u, 400u}, 10.0f);
  
  // std::size_t elems_read = nullsup::desuppress_zero(std::span<ser_type>(buffer), arr1_read);

  // std::cout << "read " << elems_read << " elements from input buffer of length " << buflen << std::endl;

  // =============
  
  // std::size_t sup_start = invec.size();
  // uint32_t num_sup = 0;

  // std::size_t ind = 0;
  // std::size_t ind_out = 0;
  // while(ind < invec.size()) {

  //   if(invec[ind] == to_suppress) {
  //     sup_start = ind;
  //     num_sup = 0;
      
  //     while(invec[++ind] == to_suppress) {
  // 	num_sup++;
  //     }

  //     buffer[ind_out++] = reinterpret_cast<const ser_type&>(to_suppress);
  //     buffer[ind_out++] = num_sup;      
  //   }
  //   else {
  //     buffer[ind_out++] = reinterpret_cast<const ser_type&>(invec[ind++]);      
  //   }
  // }
  
  // std::filesystem::path testpath = "./testVector.bin";
  
  // std::fstream ofs;
  // ofs.open(testpath, std::ios::out | std::ios::binary);  
  // stor::Traits<std::vector<float>>::serialize(ofs, invec);
  // ofs.close();

  // std::fstream ifs;
  // ifs.open(testpath, std::ios::in | std::ios::binary);
  // std::vector<float> vec_read = stor::Traits<std::vector<float>>::deserialize(ifs);
  // ifs.close();

  // for(std::size_t ind = 0; ind < invec.size(); ind++) {
  //   if(invec[ind] != vec_read[ind]) {
  //     std::cout << "problem!" << std::endl;
  //   }
  // }

  std::cout << "done" << std::endl;
}
