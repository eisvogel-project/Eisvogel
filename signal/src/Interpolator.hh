#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include <cmath>
#include <iostream>

#include "Common.hh"
#include "Kernels.hh"
#include "NDArray.hh"
#include "IteratorUtils.hh"

template <template<typename, std::size_t> class ArrayT, typename ValueT, std::size_t dims>
class Interpolator {

public:

  Interpolator(const ArrayT<ValueT, dims>& data, const Kernel& kernel) : 
    m_data(data), m_kernel(kernel) { };

  template <typename... FracInds>
  ValueT Interpolate(FracInds... frac_inds) const requires(sizeof...(FracInds) == dims) {
    DenseVector<scalar_t> target_inds({static_cast<scalar_t>(frac_inds)...});
    return Interpolate(target_inds);
  }

  ValueT Interpolate(const DenseVector<scalar_t>& target_inds) const {
    IndexVector start_inds(dims, 0);
    IndexVector end_inds(dims, 0);

    for(std::size_t i = 0; i < dims; i++) {
      scalar_t start_ind = std::ceil(target_inds(i) - m_kernel.Support());
      scalar_t end_ind = std::floor(target_inds(i) + m_kernel.Support() + 1);

      // extrapolation is not permitted
      if((start_ind < 0) || (end_ind > m_data.shape(i))) {
	std::cerr << "Extrapolation is not permitted" << std::endl;
	throw;
      }

      start_inds(i) = std::size_t(start_ind);
      end_inds(i) = std::size_t(end_ind);
    }

    ValueT interpolated_value = ValueT();

    // std::cout << "--------------" << std::endl;
    // std::cout << "interpolating onto: ";
    // for(auto cur: target_inds) {
    //   std::cout << cur << "  ";
    // }
    // std::cout << std::endl;

    // iterate over all dimensions
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      ValueT kernel_weight = 1.0;
      for(std::size_t i = 0; i < dims; i++) {
	kernel_weight *= m_kernel(target_inds(i) - cnt(i));
      }
      interpolated_value += m_data(cnt.index()) * kernel_weight;      

      // std::cout << "cur pt = ";
      // for(auto cur: cnt.index()) {
      // 	std::cout << cur << "  ";
      // }
      // std::cout << " --> val = " << m_data(cnt.index()) << ", weight = " << kernel_weight << std::endl;
    }

    // std::cout << "interpolated value: " << interpolated_value << std::endl;
    // std::cout << "--------------" << std::endl;
    
    return interpolated_value;
  }

private:

  const ArrayT<ValueT, dims>& m_data;
  const Kernel& m_kernel;  
};

#endif
