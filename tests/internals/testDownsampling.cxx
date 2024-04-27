#include <iostream>
#include <cassert>
#include <random>

#include "Vector.hh"
#include "NDVecArray.hh"
#include "NDVecArrayOperations.hh"

template <typename T, std::size_t dims, std::size_t vec_dims, class CallableT>
void fill_array(NDVecArray<T, dims, vec_dims>& to_fill, const Vector<std::size_t, dims>& start_ind, CallableT&& fill_function) {

  auto worker = [&](const Vector<std::size_t, dims>& ind) {
    to_fill[ind] = fill_function(start_ind + ind);
  };  
  IteratorUtils::index_loop_over_array_elements(to_fill, worker);
}

template <typename T, std::size_t dims, std::size_t vec_dims>
Vector<T, vec_dims> ind_sum(const Vector<std::size_t, dims>& ind) {

  std::mt19937 gen64;
  
  std::size_t hash = 0;
  std::size_t mulfact = 1;
  for(std::size_t i = 0; i < dims; i++) {
    hash += ind[i] * mulfact;
    mulfact *= 100;
  }
  gen64.seed(hash);
  
  std::size_t sum = 0;
  for(std::size_t i = 0; i < dims; i++) {
    sum += ind[i] * i;
  }

  Vector<T, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    float randval = static_cast <float> (gen64()) / static_cast <float> (gen64.max());
    retvec[i] = std::sin(sum + i) + randval;
  }
  
  return retvec;
}

template <std::size_t dims, std::size_t vec_dims>
void run_test(const Vector<std::size_t, dims>& original_shape, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& downsampling_factor) {

  auto filler = [&](const Vector<std::size_t, dims>& ind) {
    return ind_sum<float, dims, vec_dims>(ind);
  };
  
  NDVecArray<float, dims, vec_dims> to_downsample(original_shape);
  fill_array(to_downsample, start_ind, filler);

  Vector<std::size_t, dims> downsampled_init_shape(1);
  NDVecArray<float, dims, vec_dims> downsampled(downsampled_init_shape, 0.0f);

  Vector<std::size_t, dims> downsampled_start_ind{0u, 0u, 0u};
  Downsampling::downsample(to_downsample, start_ind, downsampling_factor, downsampled, downsampled_start_ind);

  if(downsampled.GetNumberElements() == 0) {
    std::cout << "Got empty downsampled array" << std::endl;    
  }
  
  Vector<std::size_t, dims> expected_downsampled_shape = Downsampling::get_downsampled_shape(start_ind, original_shape, downsampling_factor);
  if(expected_downsampled_shape != downsampled.GetShape()) {
    std::cout << "Problem: expected downsampled shape = " << expected_downsampled_shape << ", got " << downsampled.GetShape() << std::endl;
    throw;
  }
  
  std::cout << "Got downsampled array with shape = " << downsampled.GetShape() << ", with downsampled_start_ind = " << downsampled_start_ind << std::endl;

  auto checker = [&](const Vector<std::size_t, dims>& downsampled_ind) {

    Vector<std::size_t, dims> global_downsampled_ind = downsampled_ind + downsampled_start_ind;
    Vector<std::size_t, dims> global_ind = global_downsampled_ind * downsampling_factor;

    assert(to_downsample.has_index(global_ind - start_ind));
    assert(downsampled.has_index(global_downsampled_ind - downsampled_start_ind));
    
    for(std::size_t i = 0; i < vec_dims; i++) {   
      if(to_downsample[global_ind - start_ind][i] != downsampled[global_downsampled_ind - downsampled_start_ind][i]) {
	std::cout << "Problem at original = " << global_ind << ", downsampled = " << global_downsampled_ind << std::endl;
	throw;
      }
    }
  };
  IteratorUtils::index_loop_over_array_elements(downsampled, checker);
}

int main(void) {

  constexpr std::size_t dims = 3;
  constexpr std::size_t vec_dims = 2;
  
  Vector<std::size_t, dims> downsampling_factor{1u, 2u, 4u};
  Vector<std::size_t, dims> original_shape{1u, 170u, 20u};
  Vector<std::size_t, dims> start_ind{0u, 1u, 4u};  

  run_test<dims, vec_dims>(original_shape, start_ind, downsampling_factor);
  
  std::cout << "done" << std::endl;
}
