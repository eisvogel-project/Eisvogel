#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>

#include "NDVecArray.hh"
#include "DistributedNDVecArray.hh"
#include "Eisvogel/IteratorUtils.hh"

using T = float;

template <std::size_t dims, std::size_t vec_dims>
Vector<T, vec_dims> ind_sum(const Vector<std::size_t, dims>& ind) {

  std::size_t sum = 0;
  for(const std::size_t& cur : ind) {
    sum += cur;
  }

  Vector<T, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    retvec[i] = std::sin(sum + i);
  }
  
  return retvec;
}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void fill_array(NDVecArray<T, dims, vec_dims>& to_fill, const Vector<std::size_t, dims>& start_ind, CallableT&& fill_function) {

  auto worker = [&](const Vector<std::size_t, dims>& ind) {
    to_fill[ind] = fill_function(start_ind + ind);
  };  
  index_loop_over_array_elements(to_fill, worker);
}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void register_chunks(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
		     const Vector<std::size_t, dims>& chunk_size, CallableT&& filler) {

  auto chunk_registerer = [&](const Vector<std::size_t, dims>& chunk_start, const Vector<std::size_t, dims>& chunk_end) {

    // prepare filled array
    NDVecArray<T, dims, vec_dims> array_buffer(chunk_size);    
    fill_array(array_buffer, chunk_start, filler);

    // register as chunk
    std::cout << "Registering chunk: " << chunk_start << " -> " << chunk_end << std::endl;
    darr.RegisterChunk(chunk_start, array_buffer);
  };
  
  index_loop_over_chunks(start_ind, end_ind, chunk_size, chunk_registerer);  
}

template <std::size_t dims, std::size_t vec_dims>
void append_slices(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, std::size_t number_slices_to_append) {

}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void test_correctness(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, CallableT&& filler) {

  using view_t = typename DistributedNDVecArray<NDVecArray, T, dims, vec_dims>::view_t;
  auto checker = [&](const Vector<std::size_t, dims>& ind) {
    
    view_t darr_elems = darr[ind];
    Vector<T, vec_dims> filled_values = filler(ind);
    
    for(std::size_t i = 0; i < vec_dims; i++) {

      T darr_value = darr_elems[i];
      T filled_value = filled_values[i];
      
      if(darr_value != filled_value) {
	std::cout << "Error: discrepancy for index " << ind << ", found " << darr_value << ", expected " << filled_value << "!" << std::endl;
	throw;
      }
    }
  };

  Vector<std::size_t, dims> start_ind(0);
  Vector<std::size_t, dims> shape = darr.GetShape();  
  index_loop_over_elements(start_ind, shape, checker);
}

int main(int argc, char* argv[]) {

  constexpr std::size_t dims = 3;
  constexpr std::size_t vec_dims = 2;
  using darr_t = DistributedNDVecArray<NDVecArray, T, dims, vec_dims>;
    
  std::filesystem::path workdir = "./darr_test";
  if(std::filesystem::exists(workdir)) {
    std::filesystem::remove_all(workdir);
  }

  Vector<std::size_t, dims> init_cache_el_shape(100);
  Vector<std::size_t, dims> streamer_chunk_shape(stor::INFTY);
  streamer_chunk_shape[0] = 1;
  
  darr_t darr(workdir, 2, init_cache_el_shape, streamer_chunk_shape);
  
  Vector<std::size_t, dims> chunk_size(400);
  Vector<std::size_t, dims> start_ind(0);
  Vector<std::size_t, dims> end_ind(1000);

  auto filler = [&](const Vector<std::size_t, dims>& ind){return ind_sum<dims, vec_dims>(ind);};
  
  register_chunks(darr, start_ind, end_ind, chunk_size, filler);
  test_correctness(darr, filler);
  
  std::cout << "done" << std::endl;
}
