#include <iostream>
#include <fstream>
#include <vector>

#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"
#include "SignalExport.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string gf_path = argv[1];
  SignalCalculator calc(gf_path, 50);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t z_ant = -580;
  scalar_t bz = 300;
  scalar_t bx = 250;
  scalar_t z0_ice = 10;
  scalar_t tstart = -250, tend = 250;
  scalar_t beta = 0.99;
  auto charge = [&](scalar_t t){
    scalar_t cur_x = bx - beta * t;
    return -std::exp(-std::pow((cur_x - bx) / z0_ice, 2.0));
  };

  scalar_t t_sig_start = 650;
  scalar_t t_sig_samp = 0.5;
  std::size_t num_samples = 400;
  std::vector<scalar_t> signal_values(num_samples);
  std::vector<scalar_t> signal_times(num_samples);

  for(std::size_t ind = 0; ind < num_samples; ind++) {
    signal_times[ind] = t_sig_start + ind * t_sig_samp;
  }

  std::cout << "Integrating signal ..." << std::endl;

  std::vector<scalar_t> x_values;
  std::vector<scalar_t> ce_values;
  
  scalar_t delta_t = 1.0f;
  for(scalar_t cur_t = tstart; cur_t < tend; cur_t += delta_t) {
    
    scalar_t cur_x = bx - beta * cur_t;
    LineCurrentSegment cur_track(XYZCoordVector{bx - beta * cur_t, 0.0f, z_ant + bz},   // track start position
				 XYZCoordVector{bx - beta * (cur_t + delta_t), 0.0f, z_ant + bz},   // track end position
				 cur_t, cur_t + delta_t, charge(cur_t));    

    x_values.push_back(cur_x);
    ce_values.push_back(charge(cur_t));
    
    // std::cout << "charge at t = " << cur_t << ": " << charge(cur_t) << std::endl;
    
    calc.AccumulateSignal(cur_track, t_sig_start, t_sig_samp, num_samples, signal_values);       
  }

  ExportSignal(signal_times, signal_values, "signal.csv");
  ExportSignal(x_values, ce_values, "ce.csv");
  
  for(scalar_t& cur: signal_values) {
    std::cout << cur << std::endl;
  }
  
  return 0;
}
