#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include <cmath>

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
      start_inds(i) = std::size_t(std::clamp(std::ceil(target_inds(i) - m_kernel.Support()), 
					     scalar_t(0), scalar_t(m_data.shape(i)))
				  );
      end_inds(i) = std::size_t(std::clamp(std::floor(target_inds(i) + m_kernel.Support() + 1), 
					   scalar_t(0), scalar_t(m_data.shape(i)))
				);
    }

    ValueT interpolated_value = ValueT();
    IndexVector cur_inds = start_inds;

    // iterate over all dimensions
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      ValueT kernel_weight = 1.0;
      for(std::size_t i = 0; i < dims; i++) {
	kernel_weight *= m_kernel(target_inds(i) - cnt(i));
      }
      interpolated_value += m_data(cnt.index()) * kernel_weight;      
    }   
    
    return interpolated_value;
  }

private:

  const ArrayT<ValueT, dims>& m_data;
  const Kernel& m_kernel;  
};

#endif
