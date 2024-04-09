namespace Integrator {

  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, const Current0D& curr, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal) {
    
    // 1) order current segments
    
    // 2) calculate and accumulate the signal for each segment
    
  }
  
  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, std::vector<LineCurrentSegment>& curr_segs, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal) {

    // some logic here to sort the `curr_segs` so that spatially neighbouring segments are evaluated successively and long segments are broken up
    
  }
  
  template <class GreensFunctionT>
  void accumulate_integrate(GreensFunction& gf, const LineCurrentSegment& curr_seg, scalar_t t_start, scalar_t t_end, scalar_t t_samp, std::vector<scalar_t>& signal) {

    // within the context of a single line segment, don't need to be smart but integrate right away
    
  }
}
