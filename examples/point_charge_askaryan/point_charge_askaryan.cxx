#include <iostream>
#include <fstream>
#include <vector>

#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Current0DOld.hh"
#include "Eisvogel/SignalExport.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];
  SignalCalculator calc(wf_path);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t z0 = 100;
  scalar_t tstart = -250, tend = 250;
  scalar_t beta = 0.9;
  auto charge = [&](scalar_t t){
    return -exp(beta * t / z0);
  };

  std::cout << "Building trajectory ..." << std::endl;

  std::vector<CoordVector> points;
  std::vector<scalar_t> charges;

  for(scalar_t cur_t = tstart; cur_t < tend; cur_t += 1) {
    points.push_back(CoordUtils::MakeCoordVectorTXYZ(cur_t, beta * cur_t, 0, b));
    charges.push_back(charge(cur_t));
  }
  points.push_back(CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b));
  Current0D track(std::move(points), std::move(charges));

  std::cout << "Computing signal ..." << std::endl;
  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = -10; cur_t < 20; cur_t += 1) {
    scalar_t cur_signal = calc.ComputeSignal(track, cur_t);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
