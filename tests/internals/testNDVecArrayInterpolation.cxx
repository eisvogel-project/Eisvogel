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

// evaluates a linear function
template <std::size_t dims, std::size_t vec_dims>
Vector<scalar_t, vec_dims> linear(const Vector<scalar_t, dims>& pos, const Vector<scalar_t, dims>& coeffs) {

  scalar_t lin_val = VectorUtils::inner(pos, coeffs);  
  Vector<scalar_t, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    retvec[i] = (i + 1) * lin_val;
  }
  
  return retvec;
}

// Fill the array with values from a linear function
template <std::size_t dims, std::size_t vec_dims>
void fill_array_linear(NDVecArray<T, dims, vec_dims>& to_fill, const Vector<scalar_t, dims>& coeffs) {

  auto linear_filler = [&](const Vector<std::size_t, dims>& ind) {    
    return linear<dims, vec_dims>(ind.template as_type<scalar_t>(), coeffs);
  };

  fill_array(to_fill, linear_filler);  
}

template <class KernelT, std::size_t dims, std::size_t vec_dims>
void test_closure_linear_data(const Vector<std::size_t, dims>& data_shape, const Vector<scalar_t, dims>& coeffs, Vector<scalar_t, dims-1> outer_inds,
			      T inner_ind_start, T inner_ind_end, T inner_ind_delta, scalar_t rel_th = 1e-6) {

  NDVecArray<T, dims, vec_dims> data(data_shape);
  fill_array_linear<dims, vec_dims>(data, coeffs);
  
  // prepare space for the interpolated result
  std::size_t num_inner_inds = (std::size_t)((inner_ind_end - inner_ind_start) / inner_ind_delta) + 1;
  
  Vector<std::size_t, 1> shape_result(num_inner_inds);
  NDVecArray<T, 1, vec_dims> interp_vals(shape_result, 0.0);

  // perform interpolation
  auto start = std::chrono::high_resolution_clock::now();
  Interpolation::interpolate<KernelT>(data, interp_vals, outer_inds, inner_ind_start, inner_ind_end, inner_ind_delta);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);  

  // test closure
  Vector<scalar_t, dims> cur_pos;
  for(std::size_t i = 0; i < dims-1; i++) {
    cur_pos[i] = outer_inds[i];
  }

  std::size_t elem_ind = 0;
  for(T cur_inner = inner_ind_start; cur_inner < inner_ind_end; cur_inner += inner_ind_delta) {
    cur_pos[dims-1] = cur_inner;
    Vector<scalar_t, vec_dims> expected_value = linear<dims, vec_dims>(cur_pos, coeffs);
    auto interpolated_value = interp_vals[elem_ind];

    for(std::size_t i = 0; i < vec_dims; i++) {
      if(std::fabs((interpolated_value[i] - expected_value[i]) / expected_value[i]) > rel_th) {
	std::cout << "Found non-closure!" << std::endl;
	throw;
      }
    }    
    elem_ind++;
  }
  std::cout << "Performed " << elem_ind << " interpolations in " << duration << " (" << 1e6 * elem_ind / duration.count() << " interpolations/sec); all OK" << std::endl;
}

int main(int argc, char* argv[]) {
  
  // fill array with linear data
  Vector<std::size_t, dims> shape(400u);
  Vector<scalar_t, dims> coeffs(1.0f);

  // interpolation target
  Vector<T, dims-1> outer_inds{250.123f, 250.567f};  
  T inner_ind_start = 2.123;
  T inner_ind_end = 350.15;
  T inner_ind_delta = 0.5123;

  test_closure_linear_data<Interpolation::Kernel::Keys, dims, vec_dims>(shape, coeffs, outer_inds, inner_ind_start, inner_ind_end, inner_ind_delta);

  test_closure_linear_data<Interpolation::Kernel::Linear, dims, vec_dims>(shape, coeffs, outer_inds, inner_ind_start, inner_ind_end, inner_ind_delta);
  
  std::cout << "done" << std::endl;
}
