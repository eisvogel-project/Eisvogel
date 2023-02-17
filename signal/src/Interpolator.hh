#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include "Common.hh"
#include "NDArray.hh"
#include "IteratorUtils.hh"

template <template<typename, std::size_t> class ArrayT, typename ValueT, std::size_t dims>
class Interpolator {

public:

  Interpolator(const ArrayT<ValueT, dims>& data) : m_data(data) { };

  template <typename... FracInds>
  ValueT Interpolate(FracInds... frac_inds) requires(sizeof...(FracInds) == dims) {

    std::size_t support = 5;

    DenseVector<scalar_t> target_inds({static_cast<scalar_t>(frac_inds)...});
    IndexVector start_inds(target_inds.size(), 0);
    IndexVector end_inds(target_inds.size(), 0);

    for(std::size_t i = 0; i < target_inds.size(); i++) {
      start_inds(i) = std::size_t(std::clamp(target_inds(i) - support, scalar_t(0), scalar_t(m_data.shape(i))));
      end_inds(i) = std::size_t(std::clamp(target_inds(i) + support, scalar_t(0), scalar_t(m_data.shape(i))));
    }

    ValueT interpolated_value = ValueT();
    IndexVector cur_inds = start_inds;

    std::cout << "------ begin iter" << std::endl;

    // iterate over all dimensions
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
      for(size_t cur: cnt.index()) {
	std::cout << cur << "  ";
      }
      std::cout << std::endl;
    }
    
    std::cout << "------ end iter" << std::endl;   
    
    return interpolated_value;
  }

private:

  const ArrayT<ValueT, dims>& m_data;

};

#endif
