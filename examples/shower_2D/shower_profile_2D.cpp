#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/SignalExport.hh"
#include "shower_creator.h"
#include "shower_2D.h"
#include "constants.h"
#include "units.h"

namespace CU = CoordUtils;

int main(int argc, char* argv[]) {

  std::string wf_path = argv[1];
  SignalCalculator signal_calc(wf_path);
			

  std::array<float, 3> shower_vertex = {-346 * 2, .1, -256 * 2};
  std::cout << "Reading Showers \n";

  showers::ShowerCreator2D shower_creator("/home/welling/RadioNeutrino/scripts/Eisvogel/extern/shower_profile/shower_file");
    std::cout << "Building shower \n";
  showers::Shower2D new_shower = shower_creator.create_shower(
          shower_vertex,
          5.0e+18,
          90 * units::degree,
          0,
          0
          );
  new_shower.dump_profile();
}
