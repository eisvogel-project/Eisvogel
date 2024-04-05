#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "NDVecArray.hh"
#include "DistributedNDVecArray.hh"
#include "Eisvogel/IteratorUtils.hh"

using T = float;

template <std::size_t dims, std::size_t vec_dims>
Vector<T, vec_dims> ind_sum(const Vector<std::size_t, dims>& ind) {

  std::size_t hash = 0;
  std::size_t mulfact = 1;
  for(std::size_t i = 0; i < dims; i++) {
    hash += ind[i] * mulfact;
    mulfact *= 100;
  }
  srand(hash);
  
  std::size_t sum = 0;
  for(std::size_t i = 0; i < dims; i++) {
    sum += ind[i] * i;
  }

  Vector<T, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    float randval = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    retvec[i] = std::sin(sum + i) + randval;
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
		     const Vector<std::size_t, dims>& chunk_shape, CallableT&& filler) {

  auto chunk_registerer = [&](const Vector<std::size_t, dims>& chunk_start, const Vector<std::size_t, dims>& chunk_end) {

    // prepare filled array
    NDVecArray<T, dims, vec_dims> array_buffer(chunk_shape);    
    fill_array(array_buffer, chunk_start, filler);

    // register as chunk
    std::cout << "Registering chunk: " << chunk_start << " -> " << chunk_end << std::endl;
    darr.RegisterChunk(array_buffer, chunk_start);
  };
  
  index_loop_over_chunks(start_ind, end_ind, chunk_shape, chunk_registerer);  
}

template <std::size_t axis, std::size_t dims, std::size_t vec_dims, class CallableT>
void append_slices(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, const Vector<std::size_t, dims>& slice_shape, CallableT&& filler) {

  Vector<std::size_t, dims> shape = darr.GetShape();

  // start index for slice concatenation
  Vector<std::size_t, dims> start_ind(0);
  start_ind[axis] = shape[axis];

  // end index of the extended array
  Vector<std::size_t, dims> end_ind = shape;
  end_ind[axis] += slice_shape[axis];
  
  auto chunk_appender = [&](const Vector<std::size_t, dims>& chunk_start, const Vector<std::size_t, dims>& chunk_end) {

    // prepare filled array
    NDVecArray<T, dims, vec_dims> array_buffer(slice_shape);
    fill_array(array_buffer, chunk_start, filler);

    // append as chunk
    std::cout << "Appending slice at " << chunk_start << std::endl;
    darr.template AppendSlice<axis>(chunk_start, array_buffer);
  };
  
  index_loop_over_chunks(start_ind, end_ind, slice_shape, chunk_appender);
}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void test_darr_correctness(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, CallableT&& filler) {

  std::cout << "testing closure ... ";
  
  using view_t = typename DistributedNDVecArray<NDVecArray, T, dims, vec_dims>::view_t;
  auto checker = [&](const Vector<std::size_t, dims>& ind, const view_t& darr_elems) {
    
    Vector<T, vec_dims> filled_values = filler(ind);
    
    for(std::size_t i = 0; i < vec_dims; i++) {

      T darr_value = darr_elems[i];
      T filled_value = filled_values[i];
      
      if(darr_value != filled_value) {
	std::cout << "Error: discrepancy for element index " << ind << ", vector index " << i << ": found " << darr_value
		  << ", expected " << filled_value << "!" << std::endl;
	throw;
      }
    }
  };
  darr.index_loop_over_elements(checker);

  std::cout << "OK!" << std::endl;
}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void test_fill_array(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr,
		     const Vector<std::size_t, dims>& region_start_ind, const Vector<std::size_t, dims>& region_end_ind,
		     CallableT&& filler) {

  std::cout << "testing array filling ... ";
  
  // construct a buffer with the neede shape
  NDVecArray<T, dims, vec_dims> extracted_array(region_end_ind - region_start_ind);
  
  darr.FillArray(extracted_array, region_start_ind, region_end_ind);

  using view_t = typename NDVecArray<T, dims, vec_dims>::view_t;
  auto checker = [&](const Vector<std::size_t, dims>& ind) {

    Vector<T, vec_dims> filled_values = filler(ind);
    view_t extracted_values = extracted_array[ind - region_start_ind];
    
    for(std::size_t i = 0; i < vec_dims; i++) {

      T filled_value = filled_values[i];
      T extracted_value = extracted_values[i];

      if(filled_value != extracted_value) {
	std::cout << "Error: discrepancy for element index " << ind << ", vector index " << i << ": found " << extracted_value
		  << ", expected " << filled_value << "!" << std::endl;
	throw;
      }      
    }
  };
  
  index_loop_over_elements(region_start_ind, region_end_ind, checker);

  std::cout << "OK!" << std::endl;
}

int main(int argc, char* argv[]) {

  constexpr std::size_t dims = 3;
  constexpr std::size_t vec_dims = 2;
  using darr_t = DistributedNDVecArray<NDVecArray, T, dims, vec_dims>;
    
  std::filesystem::path workdir = "./darr_test";
  if(std::filesystem::exists(workdir)) {
    std::filesystem::remove_all(workdir);
  }

  Vector<std::size_t, dims> init_cache_el_shape(10);
  Vector<std::size_t, dims> streamer_chunk_shape(stor::INFTY);
  streamer_chunk_shape[0] = 1;
  
  darr_t darr(workdir, 2, init_cache_el_shape, streamer_chunk_shape);
  
  Vector<std::size_t, dims> chunk_size(3);
  Vector<std::size_t, dims> start_ind(0);
  Vector<std::size_t, dims> end_ind(11);

  auto filler = [&](const Vector<std::size_t, dims>& ind){return ind_sum<dims, vec_dims>(ind);};
  
  register_chunks(darr, start_ind, end_ind, chunk_size, filler);

  Vector<std::size_t, dims> slice_shape(3);
  slice_shape[0] = 10;
  
  append_slices<0>(darr, slice_shape, filler);
  append_slices<0>(darr, slice_shape, filler);
  
  test_darr_correctness(darr, filler);

  std::cout << darr.GetShape() << std::endl;
  
  darr.SwapAxes<0, 2>();

  std::cout << darr.GetShape() << std::endl;

  auto swapped_filler = [&](const Vector<std::size_t, dims>& ind){
    Vector<std::size_t, dims> swapped_ind = ind;
    std::swap(swapped_ind[0], swapped_ind[2]);    
    return ind_sum<dims, vec_dims>(swapped_ind);
  };  

  test_darr_correctness(darr, swapped_filler);

  auto boundary_evaluator = [](){};
  
  std::filesystem::path workdir_tmp = "./darr_test_tmp";
  Vector<std::size_t, dims> requested_chunk_size(4);  
  darr.RebuildChunks(requested_chunk_size, workdir_tmp, 1, boundary_evaluator);

  std::cout << darr.GetShape() << std::endl;
  
  test_darr_correctness(darr, swapped_filler);
    
  std::cout << "done" << std::endl;
}
