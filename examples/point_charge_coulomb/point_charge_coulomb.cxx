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

  std::string wf_path = argv[1];
  SignalCalculator calc(wf_path);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = -250, tend = 250;
  scalar_t charge = 1;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;  
  LineCurrentSegment track(XYZCoordVector{beta * tstart, 0.0f, b},    // track start position
			   XYZCoordVector{beta * tend, 0.0f, b},   // track end position
			   tstart, tend, charge);

  scalar_t t_sig_start = -10;
  scalar_t t_sig_samp = 1.0;
  std::size_t num_samples = 30;
  std::vector<scalar_t> signal_values(num_samples);
  
  std::cout << "Computing signal ..." << std::endl;

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
