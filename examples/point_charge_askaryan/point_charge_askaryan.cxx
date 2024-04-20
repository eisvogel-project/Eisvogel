#include <iostream>
#include <fstream>
#include <vector>

#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string gf_path = argv[1];
  SignalCalculator calc(gf_path);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t z0 = 100;
  scalar_t tstart = -250, tend = 250;
  scalar_t beta = 0.9;
  auto charge = [&](scalar_t t){
    return -exp(beta * t / z0);
  };

  scalar_t t_sig_start = -10;
  scalar_t t_sig_samp = 1.0;
  std::size_t num_samples = 30;
  std::vector<scalar_t> signal_values(num_samples);

  std::cout << "Integrating signal ..." << std::endl;
  
  scalar_t delta_t = 1.0f;
  for(scalar_t cur_t = tstart; cur_t < tend; cur_t += delta_t) {

    LineCurrentSegment cur_track(XYZCoordVector{beta * cur_t, 0.0f, b},               // track start position
				 XYZCoordVector{beta * (cur_t + delta_t), 0.0f, b},   // track end position
				 cur_t, cur_t + delta_t, charge(cur_t));
 
    calc.AccumulateSignal(cur_track, t_sig_start, t_sig_samp, num_samples, signal_values);       
  }
  
  for(scalar_t& cur: signal_values) {
    std::cout << cur << std::endl;
  }
  
  return 0;
}
