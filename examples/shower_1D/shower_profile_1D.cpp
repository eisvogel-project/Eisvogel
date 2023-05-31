#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/SignalExport.hh"
#include "shower_creator.h"
#include "shower_1D.h"
#include "constants.h"
#include "units.h"

namespace CU = CoordUtils;

int main(int argc, char* argv[]) {

  std::string wf_path = argv[1];
  SignalCalculator signal_calc(wf_path);
			

  std::array<float, 3> shower_vertex = {-346 * 2, .1, -256 * 2};
  std::cout << "Building Shower \n";

  showers::ShowerCreator shower_creator("/home/welling/RadioNeutrino/scripts/Eisvogel/extern/shower_profile/shower_file");
  showers::Shower1D shower = shower_creator.create_shower(
          shower_vertex,
          5.0e+18,
          90 * units::degree,
          0,
          0
          );

  // Show dimensions of the required weighting field
  std::vector<double> t;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> ce;
  shower.get_shower(1., &t, &x, &y, &z, &ce);
  double r_max = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double r_min = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double z_max = z[0];
  double z_min = z[0];
  double r;
  for (int i=0; i<t.size(); i++) {
      r = sqrt(x[i] * x[i] + y[i] * y[i]);
      r_min = std::min(r_min, r);
      r_max = std::max(r_max, r);
      z_min = std::min(z_min, z[i]);
      z_max = std::max(z_max, z[i]);
  }
  r_max /= constants::c;
  r_min /= constants::c;
  z_max /= constants::c;
  z_min /= constants::c;

  std::cout << "Required weighting field size:\n";
  std::cout << r_min << "< r < " << r_max << "\n";
  std::cout << z_min << "< z < " << z_max << "\n";


  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = 0, tend = 20;
  scalar_t charge = 1;
  scalar_t beta = 0.9;
  
  
  std::vector<scalar_t> signal_times, signal_values;
  std::cout << "Integrating \n";
  Current0D current = shower.get_current(0.1);
  for(scalar_t cur_t = 3400; cur_t < 4000; cur_t += 1) {
    scalar_t cur_signal = signal_calc.ComputeSignal(current, cur_t);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }
  // convert to normal units
  for (int i; i < signal_values.size(); i++) {
      signal_values[i] = signal_values[i] / 2.218e10 * constants::c;
  }
  ExportSignal(signal_times, signal_values, "./5e18EM.csv");
  return 0;
}
