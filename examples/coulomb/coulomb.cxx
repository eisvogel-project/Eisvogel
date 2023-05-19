#include <iostream>

#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "Eisvogel/SignalExport.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string wf_path = argv[1];
  SignalCalculator calc(wf_path);

  SparseCurrentDensity3D curr_dist(CoordUtils::MakeCoordVectorTXYZ(0.1, 0.1, 0.1, 0.1));

  curr_dist.addCurrentElement({CoordUtils::MakeCoordVectorTXYZ(0, 0, 0, 0), CoordUtils::MakeFieldVectorXYZ(0, 0, 0)});

  calc.ComputeSignal(curr_dist, 0.0);

}
