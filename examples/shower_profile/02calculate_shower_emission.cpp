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

    std::string wf_path;
    std::string output_path;
    if (argc < 2) {
        wf_path = "weighting_field";
    } else {
        wf_path = argv[1];
    }
    if (argc < 3) {
        output_path = "shower_radio_emission.csv";
    } else {
        output_path = argv[2];
    }
    SignalCalculator signal_calc(wf_path);
                

    std::array<float, 3> shower_vertex = {-112  , .1, -165};
    
    showers::ShowerCreator shower_creator("/home/welling/RadioNeutrino/scripts/Eisvogel/extern/shower_profile/shower_file");
    showers::Shower1D shower = shower_creator.create_shower(
            shower_vertex,
            5.0e+18,
            90 * units::degree,
            0,
            0
            );

    // Show dimensions of the required weighting field
     shower.print_dimenstions();

    // test trajectory: a point charge moving parallel to the x-axis 
    // with a constant impact parameter of 'b' along the z-axis

    scalar_t sampling_rate = 2.;

    std::vector<scalar_t> signal_times, signal_values;
    Current0D current = shower.get_current(0.1);
    for(scalar_t cur_t = 1050; cur_t < 1300; cur_t += 1. / sampling_rate) {
        scalar_t cur_signal = signal_calc.ComputeSignal(current, cur_t);
        signal_times.push_back(cur_t);
        signal_values.push_back(cur_signal);
    }
    // convert to normal units
    for (int i; i < signal_values.size(); i++) {
        signal_values[i] = signal_values[i] / 2.218e10 * constants::c;
    }
    ExportSignal(signal_times, signal_values, output_path);
    return 0;
}
