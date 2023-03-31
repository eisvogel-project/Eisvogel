#include <iostream>
#include <fstream>
#include "Common.hh"
#include "Kernels.hh"
#include "Interpolator.hh"
#include "MathUtils.hh"
#include "CoordUtils.hh"
#include "IteratorUtils.hh"
#include "NDArray.hh"

int main(void) {

  // 1d interpolation
  DeltaVector step({0.3});
  CoordVector start_coords({1e-3});
  CoordVector end_coords({45.001});

  std::array<std::size_t, 1> number_pts;
  std::size_t dim = 0;
  for(auto cur: (end_coords - start_coords) / step) {
    number_pts[dim] = std::size_t(cur);
    dim += 1;
  }

  DenseVector sampled(number_pts, 0.0);

  DenseVector step_chosen = (end_coords - start_coords) / sampled.shape();

  IndexVector start_inds({0});
  IndexVector end_inds({number_pts});

  auto func_1d = [&](scalar_t x) -> scalar_t {
    scalar_t abs_x = std::fabs(x);
    return 1.0 / std::pow(abs_x, 3) * std::sin(abs_x);
  };

  for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
    IndexVector ind = cnt.index();

    scalar_t x = start_coords(0) + step_chosen(0) * ind(0);
    sampled(ind) = func_1d(x);

    std::cout << "sampling function at " << x << ": " << sampled(ind) << std::endl;
  }

  SplineInterpolationKernelOrder3 kernel;
  Interpolator itpl(sampled, kernel);

  for(scalar_t cur_x = 1.3; cur_x < 7.0; cur_x += 0.31) {
    CoordVector frac_inds = (CoordVector({cur_x}) - start_coords) / (end_coords - start_coords) * sampled.shape();
    std::cout << "x = " << cur_x << std::endl;
    std::cout << "f_true(x) = " << func_1d(cur_x) << std::endl;
    std::cout << "f_itpl(x) = " << itpl.Interpolate(frac_inds) << std::endl;
  }

  return 0;
}
