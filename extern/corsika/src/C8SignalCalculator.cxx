#include "C8SignalCalculator.hh"
#include "GreensFunction.hh"
#include "Current.hh"

#include <iostream>

C8SignalCalculator::C8SignalCalculator(std::filesystem::path gf_path, std::size_t cache_depth) {
  
  m_gf = std::make_shared<CylindricalGreensFunction>(gf_path, cache_depth);
  std::cout << "created EV at " << m_gf.get() << std::endl;
}

C8SignalCalculator::~C8SignalCalculator() {
  std::cout << "EVEVEVEV C8SignalCalculator destructor EVEVEVEVEV" << std::endl;
}

void C8SignalCalculator::calculate(const std::array<float, 4>& track_start_xyzt, const std::array<float, 4>& track_end_xyzt, float track_charge,
				   float t_sig_start, float t_sig_samp, std::size_t num_samples, std::vector<float>& signal) {

  LineCurrentSegment track(XYZCoordVector{track_start_xyzt[0], track_start_xyzt[1], track_start_xyzt[2]},    // track start position
			   XYZCoordVector{track_end_xyzt[0], track_end_xyzt[1], track_end_xyzt[2]},   // track end position
			   track_start_xyzt[3], track_end_xyzt[3], track_charge);
    
  m_gf -> apply_accumulate<Interpolation::Kernel::Keys>(track, t_sig_start, t_sig_samp, num_samples, signal);
}
