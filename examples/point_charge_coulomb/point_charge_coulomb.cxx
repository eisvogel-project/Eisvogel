#include <iostream>
#include <fstream>

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

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = -250, tend = 250;
  scalar_t charge = 1;
  scalar_t beta = 0.9;

  std::cout << "Building trajectory ..." << std::endl;
  Current0D track({
      CoordUtils::MakeCoordVectorTXYZ(tstart, beta * tstart, 0, b),
  	CoordUtils::MakeCoordVectorTXYZ(tend, beta * tend, 0, b)
  	},
    {charge}
    );

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
