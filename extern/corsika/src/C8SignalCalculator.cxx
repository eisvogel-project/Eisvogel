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

  std::cout << "====================" << std::endl;
  std::cout << "EV track start (x, y, z, t) = " << track_start_xyzt[0] << ", " << track_start_xyzt[1] << ", " << track_start_xyzt[2] << ", " << track_start_xyzt[3] << std::endl;
  std::cout << "EV track end (x, y, z, t) = " << track_end_xyzt[0] << ", " << track_end_xyzt[1] << ", " << track_end_xyzt[2] << ", " << track_end_xyzt[3] << std::endl;
  std::cout << "EV track charge = " << track_charge << std::endl;
  std::cout << "EV t_sig_start = " << t_sig_start << ", t_sig_samp = " << t_sig_samp << ", number_samples = " << num_samples << std::endl;
  std::cout << "====================" << std::endl;

  LineCurrentSegment track(XYZCoordVector{track_start_xyzt[0], track_start_xyzt[1], track_start_xyzt[2]},    // track start position
			   XYZCoordVector{track_end_xyzt[0], track_end_xyzt[1], track_end_xyzt[2]},   // track end position
			   track_start_xyzt[3], track_end_xyzt[3], track_charge);
  
  // // -----------

  // scalar_t b = 10;
  // scalar_t tstart = -250, tend = 250;
  // scalar_t charge = 1;
  // scalar_t beta = 0.9;

  // std::cout << "Building trajectory ..." << std::endl;  
  // LineCurrentSegment track(XYZCoordVector{beta * tstart, 0.0f, b},    // track start position
  // 			   XYZCoordVector{beta * tend, 0.0f, b},   // track end position
  // 			   tstart, tend, charge);

  // t_sig_start = -10;
  // t_sig_samp = 1.0;
  // num_samples = 30;

  // signal.resize(num_samples);
  
  // // -----------

  std::cout << "EV starting" << std::endl;
  std::cout << "EV at: " << m_gf.get() << std::endl;
  std::cout << "EV calling into apply_accumulate" << std::endl;
  
  m_gf -> apply_accumulate<Interpolation::Kernel::Keys>(track, t_sig_start, t_sig_samp, num_samples, signal);
  std::cout << "EV done" << std::endl;
}
