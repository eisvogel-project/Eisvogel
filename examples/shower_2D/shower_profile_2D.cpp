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

  std::ofstream output_file;
  output_file.open("current.csv");

  std::array<float, 3> shower_vertex = {-112, .1, -165};
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
  std::cout << "Getting Current \n";  
  DeltaVector voxel_size = CU::MakeDeltaVectorTXYZ(.2, .2, .2, .2);
  SparseCurrentDensity3D current_density = new_shower.get_current(voxel_size);
  std::cout << "Current Elements: " << current_density.getNumberofElements() << "\n";
  std::vector<scalar_t> signal_times, signal_values;
  for(const auto& [cur_pos_txyz, cur_current_density_xyz] : current_density) {
    output_file << CU::getX(cur_pos_txyz) << ", " << CU::getY(cur_pos_txyz) << ", " << CU::getZ(cur_pos_txyz) << ", ";
    output_file << CU::getXComponent(cur_current_density_xyz) << ", " << CU::getYComponent(cur_current_density_xyz) << ", " << CU::getZComponent(cur_current_density_xyz) << "\n";
  }
  std::cout << "\n";
  for(scalar_t cur_t = 1050; cur_t < 1300; cur_t += .1) {
    scalar_t cur_signal = signal_calc.ComputeSignal(current_density, cur_t);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }
    for (int i; i < signal_values.size(); i++) {
      signal_values[i] = signal_values[i] / 2.218e10 * constants::c;
  }
  ExportSignal(signal_times, signal_values, "./5e18EM_2D.csv");


}
