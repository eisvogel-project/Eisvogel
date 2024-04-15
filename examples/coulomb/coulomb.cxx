#include <iostream>

#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculatorOld.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "Eisvogel/SignalExport.hh"
#include "Eisvogel/CoordUtils.hh"

namespace CU = CoordUtils;

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];
  SignalCalculatorOld calc(wf_path);

  DeltaVector voxel_size = CU::MakeCoordVectorTXYZ(0.1, 0.1, 0.1, 0.1);

  // test trajectory: a point charge moving parallel to the x-axis 
  // with a constant impact parameter of 'b' along the z-axis
  scalar_t b = 10;
  scalar_t tstart = -250, tend = 250;
  scalar_t charge = 1;
  scalar_t beta = 0.9;

  SparseCurrentDensity3D curr_dist(voxel_size);

  for(scalar_t cur_t = tstart; cur_t <= tend; cur_t += CU::getT(voxel_size)) {

    scalar_t cur_x = beta * cur_t;
    scalar_t cur_y = 0;
    scalar_t cur_z = b;
    
    scalar_t cur_jx = charge * beta / curr_dist.getVolumeElementXYZ();
    scalar_t cur_jy = 0;
    scalar_t cur_jz = 0;

    curr_dist.addCurrentElement({
	CU::MakeCoordVectorTXYZ(cur_t, cur_x, cur_y, cur_z),
	CU::MakeFieldVectorXYZ(cur_jx, cur_jy, cur_jz)
    });
  }

  std::cout << "Computing signal ..." << std::endl;
  std::vector<scalar_t> signal_times, signal_values;
  for(scalar_t cur_t = -10; cur_t < 20; cur_t += 1) {
    scalar_t cur_signal = calc.ComputeSignal(curr_dist, cur_t);
    signal_times.push_back(cur_t);
    signal_values.push_back(cur_signal);
  }

  ExportSignal(signal_times, signal_values, "./test_signal.csv");

  return 0;
}
