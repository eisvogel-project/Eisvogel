#include <iostream>
#include "Common.hh"
#include "NDArray.hh"
#include "WeightingFieldUtils.hh"
#include "Interpolator.hh"

#include "Kernels.hh"
#include "IteratorUtils.hh"

#include <cmath>
#include <tuple>

namespace WFU = WeightingFieldUtils;

// template <typename... FracInds>
// scalar_t Interpolate(FracInds... frac_inds) {

//   std::size_t support = 4;
//   std::size_t shape_start = 0;
//   std::size_t shape_end = 10;

//   std::size_t number_dims = sizeof...(FracInds);
  
//   DenseVector<scalar_t> target_inds({frac_inds...});

//   IndexVector start_inds(target_inds.size(), 0);
//   IndexVector end_inds(target_inds.size(), 0);

//   for(std::size_t i = 0; i < target_inds.size(); i++) {
//     start_inds(i) = std::size_t(std::clamp(target_inds(i) - support, scalar_t(shape_start), scalar_t(shape_end)));
//     end_inds(i) = std::size_t(std::clamp(target_inds(i) + support, scalar_t(shape_start), scalar_t(shape_end)));
//   }

//   for(auto cur: start_inds) {
//     std::cout << cur << " / ";
//   }  
//   std::cout << std::endl;

//   for(auto cur: end_inds) {
//     std::cout << cur << " / ";
//   }  
//   std::cout << std::endl;

//   IndexVector cur_inds = start_inds;

//   std::cout << "------ begin iter" << std::endl;

//   // iterate over all dimensions
//   for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
//     for(size_t cur: cnt.index()) {
//       std::cout << cur << "  ";
//     }
//     std::cout << std::endl;
//   }
  
//   std::cout << "------ end iter" << std::endl;

//   // DenseVector<std::size_t> start_inds(target_inds.size());// = target_inds - support;
//   // DenseVector<std::size_t> end_inds(target_in = target_inds + support;

//   // std::cout << central_inds(0) << std::endl;
//   // std::cout << central_inds(1) << std::endl;

//   // std::cout << start_inds(0) << std::endl;
//   // std::cout << start_inds(1) << std::endl;

//   // for(auto cur: central_inds) {
//   //   std::cout << cur << " ";
//   // }
//   // std::cout << std::endl;

//   return 0.0;
// }

int main(void) {

  unsigned int len_t = 10;

  int n = 10;

  DenseNDArray<scalar_t, 1> vector2({1, 2, 3});

  DenseWeightingField wf = WFU::CreateElectricDipoleWeightingField();

  std::cout << wf.E_r(0,0,0) << std::endl;
  std::cout << wf.E_z(0,0,0) << std::endl;
  std::cout << wf.E_phi(0,0,0) << std::endl;

  DenseNDArray<scalar_t, 3> testarr({10, 10, 10}, 2.0);
  testarr(2, 3, 1) = 18;
  IndexVector inds({2, 3, 1});

  std::cout << "HHH " << testarr(inds) << " HHH" << std::endl;

  Kernels::SplineInterpolationKernelOrder1 kernel;

  std::cout << kernel.Support() << std::endl;
  std::cout << kernel(0.78) << std::endl;
  std::cout << kernel(1.78) << std::endl;

  std::cout << "-----" << std::endl;

  DenseNDArray<scalar_t, 2> testarr2d({10, 10}, 2.0);
  Interpolator<DenseNDArray, scalar_t, 2> itpl(testarr2d);
  std::cout << itpl.Interpolate(2, 2) << std::endl;

  std::cout << "-----" << std::endl;

  return 0;
}
