#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include <cmath>
#include <iostream>
#include <type_traits>

#include "Common.hh"
#include "Kernels.hh"
#include "NDArray.hh"
#include "IteratorUtils.hh"

#include "CoordUtils.hh"

template <typename FuncT, typename KernelT,
	  typename ValueT = std::invoke_result_t<FuncT, GridVector&>>
ValueT InterpolateFunc(FuncT func, KernelT& kernel, CoordVector& target_inds) {

  std::size_t dims = target_inds.size();
  
  GridVector start_inds(dims, 0);
  GridVector end_inds(dims, 0);
  
  for(std::size_t i = 0; i < dims; i++) {
    scalar_t start_ind = std::ceil(target_inds(i) - kernel.Support());
    scalar_t end_ind = std::floor(target_inds(i) + kernel.Support() + 1);
    
    start_inds(i) = int(start_ind);
    end_inds(i) = int(end_ind);
  }

  ValueT interpolated_value = ValueT();
  
  // iterate over all dimensions
  for(GridCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
    
    scalar_t kernel_weight = 1.0;
    for(std::size_t i = 0; i < dims; i++) {
      kernel_weight *= kernel(target_inds(i) - cnt(i));
    }
    
    ValueT cur_val = func(cnt.index());
    
    interpolated_value += cur_val * kernel_weight;      
  }
  
  return interpolated_value;
}

template <template<typename, std::size_t> class ArrayT, typename ValueT, std::size_t dims>
class Interpolator {

public:

  Interpolator(ArrayT<ValueT, dims>& data, const Kernel& kernel) : 
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

    // iterate over all dimensions
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      ValueT kernel_weight = 1.0;
      for(std::size_t i = 0; i < dims; i++) {
	kernel_weight *= m_kernel(target_inds(i) - cnt(i));
      }

      scalar_t cur_val = m_data(cnt.index());
      if(std::isnan(cur_val)) {
	std::cout << "Unallowed region" << std::endl;
	throw;
      }
      interpolated_value += cur_val * kernel_weight;      
    }

    return interpolated_value;
  }

private:

  ArrayT<ValueT, dims>& m_data;
  const Kernel& m_kernel;  
};

#endif
