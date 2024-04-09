#pragma once

namespace Integrator {

  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, const Current0D& curr, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal);
  
  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, std::vector<LineCurrentSegment>& curr_segs, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal);
  
  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, const LineCurrentSegment& curr_seg, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal);  
}
  
#include "Integrator.hxx"
