#include "C8SignalCalculator.hh"
#include "GreensFunction.hh"
#include "Current.hh"
#include "Vector.hh"

#include <iostream>

C8SignalCalculator::C8SignalCalculator(std::filesystem::path gf_path, std::size_t cache_depth) {
  
  m_gf = std::make_unique<CylindricalGreensFunction>(gf_path, cache_depth);
}

C8SignalCalculator::~C8SignalCalculator() = default;

void C8SignalCalculator::get_antenna_location(std::array<float, 3>& ant_xyz) {
  // TODO: store in Green's function metadata and read back
  ant_xyz.fill(0.0f);
}

void C8SignalCalculator::calculate(const std::array<float, 4>& track_start_xyzt, const std::array<float, 4>& track_end_xyzt, float track_charge,
				   float t_sig_start, float t_sig_samp, std::size_t num_samples, std::vector<double>& signal, float weight) {
  
  LineCurrentSegment track(XYZCoordVector{track_start_xyzt[0], track_start_xyzt[1], track_start_xyzt[2]},    // track start position
			   XYZCoordVector{track_end_xyzt[0], track_end_xyzt[1], track_end_xyzt[2]},   // track end position
			   track_start_xyzt[3], track_end_xyzt[3], track_charge);

  Green::OutOfBoundsBehavior oob_mode = Green::OutOfBoundsBehavior::Ignore; 
  m_gf -> apply_accumulate<Interpolation::Kernel::Linear, double>(track, t_sig_start, t_sig_samp, num_samples, signal, oob_mode, weight);
}
