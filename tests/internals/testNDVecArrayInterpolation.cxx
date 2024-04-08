#include <iostream>
#include <chrono>
#include "Interpolation.hh"

using T = float;
constexpr std::size_t dims = 3;
constexpr std::size_t vec_dims = 2;

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void fill_array(NDVecArray<T, dims, vec_dims>& to_fill, CallableT&& fill_function) {
  auto worker = [&](const Vector<std::size_t, dims>& ind) {
    to_fill[ind] = fill_function(ind);
  };  
  IteratorUtils::index_loop_over_array_elements(to_fill, worker);
}

int main(int argc, char* argv[]) {

  // fill array with linear data
  auto linear_filler = [](const Vector<std::size_t, dims>& ind) {

    std::size_t indsum = 0;
    for(std::size_t i = 0; i < dims; i++) {
      indsum += ind[i];
    }

    Vector<T, vec_dims> retvec;
    for(std::size_t i = 0; i < vec_dims; i++) {
      retvec[i] = (i + 1) * indsum;
    }    
    return retvec;
  };

  Vector<std::size_t, dims> shape(400u);
  NDVecArray<T, dims, vec_dims> data(shape);
  fill_array(data, linear_filler);

  std::cout << data[{10u, 10u, 10u}][0] << ", " << data[{10u, 10u, 10u}][1] << std::endl;

  Vector<T, dims-1> outer_inds(0.0);
  outer_inds[0] = 250.123;
  outer_inds[1] = 250.567;

  Vector<std::size_t, 1> shape_result(5000u);
  NDVecArray<T, 1, vec_dims> interp_vals(shape_result, 0.0);
  
  T inner_ind_start = 2.123;
  T inner_ind_end = 350.15;
  T inner_ind_delta = 1.55123;

  // Interpolation::interpolate<Interpolation::Kernel::Keys>(data, interp_vals, outer_inds,
  // 							  inner_ind_start, inner_ind_end, inner_ind_delta);
  
  auto start = std::chrono::high_resolution_clock::now();
  
  Interpolation::interpolate<Interpolation::Kernel::Keys>(data, interp_vals, outer_inds,
							  inner_ind_start, inner_ind_end, inner_ind_delta);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << "took " << duration << std::endl;
  
  scalar_t outer_indsum = 0;
  for(std::size_t i = 0; i < dims - 1; i++) {
    outer_indsum += outer_inds[i];
  }
  
  std::size_t fi = 0;
  for(scalar_t f = inner_ind_start; f < inner_ind_end; f += inner_ind_delta) {
    
    std::cout << "----" << std::endl;
    std::cout << " f = " << f << std::endl;
    std::cout << interp_vals[fi][0] << ", " << interp_vals[fi][1] << std::endl;

    scalar_t indsum = outer_indsum + f;

    if(std::fabs(interp_vals[fi][0] - indsum) / indsum > 1e-4) {
      std::cout << "error" << std::endl;
      // throw;
    }

    if(std::fabs(interp_vals[fi][1] - 2 * indsum) / indsum > 1e-4) {
      std::cout << "error" << std::endl;
      // throw;
    }
    
    std::cout << indsum << ", " << 2 * indsum << std::endl;    
    std::cout << "----" << std::endl;

    fi++;
    
  }

  std::cout << "checked " << fi << " interpolations" << std::endl;

  std::cout << fi / (duration.count() * 1e-6) << " interpolations / sec" << std::endl;
  
  std::cout << "done" << std::endl;
}
