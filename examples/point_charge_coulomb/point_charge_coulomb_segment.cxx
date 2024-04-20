#include <iostream>
#include <fstream>
#include <chrono>

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
  // shows that charge conservation is automatic
  scalar_t b = 2;
  scalar_t tstart = -10, tend = 10;
  scalar_t charge = 1;
  scalar_t beta = 0.05;

  std::cout << "Building trajectory ..." << std::endl;
  LineCurrentSegment track(XYZCoordVector{b, 0.0f, beta * tstart},
			   XYZCoordVector{b, 0.0f, beta * tend},
			   tstart, tend, charge);

  scalar_t t_sig_start = -15;
  scalar_t t_sig_samp = 1.0;
  std::size_t num_samples = 45;
  
  std::cout << "Computing signal ..." << std::endl;
  std::vector<scalar_t> signal_values(num_samples);

  std::fill(signal_values.begin(), signal_values.end(), 0.0f);
  calc.AccumulateSignal(track, t_sig_start, t_sig_samp, num_samples, signal_values);
  std::fill(signal_values.begin(), signal_values.end(), 0.0f);

  auto start = std::chrono::high_resolution_clock::now();  
  calc.AccumulateSignal(track, t_sig_start, t_sig_samp, num_samples, signal_values);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << "calculation finished in " << duration << std::endl;

  for(scalar_t& cur: signal_values) {
    std::cout << cur << std::endl;
  }

  return 0;
}
