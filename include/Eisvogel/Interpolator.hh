#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include "Common.hh"
#include "IteratorUtils.hh"
#include "CoordUtils.hh"

template <typename KernelT, typename FuncT, 
	  typename ValueT = std::invoke_result_t<FuncT, GridVector&>>
ValueT InterpolateFuncNew(FuncT func, CoordVector& target_inds) {

  std::size_t dims = target_inds.size();
  
  GridVector start_inds(dims, 0);
  GridVector end_inds(dims, 0);
  
  for(std::size_t i = 0; i < dims; i++) {
    scalar_t start_ind = std::ceil(target_inds(i) - KernelT::Support);
    scalar_t end_ind = std::floor(target_inds(i) + KernelT::Support + 1);
    
    start_inds(i) = int(start_ind);
    end_inds(i) = int(end_ind);
  }

  ValueT interpolated_value = ValueT();
  
  // iterate over all dimensions
  for(GridCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
    
    scalar_t kernel_weight = 1.0;
    for(std::size_t i = 0; i < dims; i++) {
      kernel_weight *= KernelT::Weight(target_inds(i) - cnt(i));
    }
    
    ValueT cur_val = func(cnt.index());
    
    interpolated_value += cur_val * kernel_weight;      
  }
  
  return interpolated_value;
}

#endif
