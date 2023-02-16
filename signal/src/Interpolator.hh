#ifndef __INTERPOLATOR_HH
#define __INTERPOLATOR_HH

#include "NDArray.hh"

// TOOD: does caching in a sparse array make sense?
template <template<typename, std::size_t> class ArrayT, typename ValueT, std::size_t dim>
class Interpolator {

public:

  Interpolator(ArrayT<ValueT, dim>& data) : m_data(data) { };

  template <typename... FracInds>
  ValueT Interpolate(FracInds... frac_inds) {
    ValueT interpolated_value = ValueT();

    

    return interpolated_value;
  }

private:

  ArrayT<ValueT, dim>& m_data;

};

#endif
