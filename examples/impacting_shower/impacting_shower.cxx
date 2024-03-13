#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/SignalExport.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];
  SignalCalculator calc(wf_path);

  // test trajectory: a point charge moving parallel to the z-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 100;
  scalar_t zmax = 0.0;
  scalar_t zscale = 10;
  scalar_t tstart = -40, tend = 40;
  scalar_t beta = 1.0;  // for in-ice shower
  scalar_t qmax = 1e3;
  auto charge = [&](scalar_t z){
    scalar_t cur_charge = qmax * std::exp(-std::pow((z - zmax) / zscale, 2.0));
    if(cur_charge < 1) {
      cur_charge = 0;
    }
    return cur_charge;
  };

  std::cout << "Building trajectory ..." << std::endl;

  std::vector<CoordVector> points;
  std::vector<scalar_t> charges;

  for(scalar_t cur_t = tstart; cur_t < tend; cur_t += 0.1) {
    scalar_t cur_z = -beta * cur_t;
    scalar_t cur_charge = -charge(cur_z);

    std::cout << "z = " << cur_z << ", q = " << cur_charge << std::endl;
    
    points.push_back(CoordUtils::MakeCoordVectorTXYZ(cur_t, b, 0, cur_z));
    charges.push_back(cur_charge);
  }
  points.push_back(CoordUtils::MakeCoordVectorTXYZ(tend, b, 0, beta * tend));
  Current0D track(std::move(points), std::move(charges));
  
  std::cout << "Computing signal ..." << std::endl;
  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = 170; cur_t < 280; cur_t += 0.1) {
    std::cout << "calculating signal for t = " << cur_t << std::endl;
    scalar_t cur_signal = calc.ComputeSignal(track, cur_t);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./impacting_shower_signal.csv");

  return 0;
}
