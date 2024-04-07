#include <iostream>
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

  Vector<T, dims-1> outer_inds(10.123);

  Vector<std::size_t, 1> shape_result(100u);
  NDVecArray<T, 1, vec_dims> interp_vals(shape_result, 1.0);
  
  T inner_ind_start = 0.123;
  T inner_ind_end = 300.15;
  T inner_ind_delta = 10.123;

  Interpolation::interpolate<Interpolation::Kernel::Keys>(data, interp_vals, outer_inds,
			     inner_ind_start, inner_ind_end, inner_ind_delta);
  
  std::cout << "done" << std::endl;
}
