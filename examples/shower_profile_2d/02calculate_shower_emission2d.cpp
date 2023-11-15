#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "Eisvogel/Common.hh"
#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/SignalExport.hh"
#include "Eisvogel/CoordUtils.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "shower_creator.h"
#include "shower_2D.h"
#include "constants.h"
#include "units.h"

using CurrentElement = std::pair<CoordVector, FieldVector>;


int main (int argc, char* argv[]) {
    scalar_t shower_energy = 1.e18;                  // energy of the shower, in eV
    scalar_t shower_zenith = 90 * units::degree;     // zenith angle of the shower. I recommend leaving this as is
    scalar_t shower_azimuth = 0;                     // azimuth of the shower. I also recommend leaving this
    int is_hadronic = 0;                             // sets the type of the shower. Choose 1 for hardonic or 0 for electromagnetic shower
    int i_shower = 9;                                // index of the shower profile to be chosen from the library. Choose a value between 0 and 9 or set to NULL to have a shower selected at random
    std::array<float, 3> shower_vertex = {-106  , .1, -165};    // position of the shower vertex. If you change this, make sure to also adjust the weighting field accordingly
    scalar_t sampling_rate = 5.;                     // sampling rate that is used when calculating the radio signal

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
    std::cout << "Done reading weighting field \n";
    
    showers::ShowerCreator2D shower_creator("/home/welling/RadioNeutrino/scripts/corsika/plot_examples/charge_excess_profile.hdf5");
    showers::Shower2D shower = shower_creator.create_shower(shower_vertex, shower_energy, shower_zenith, shower_azimuth, is_hadronic);
    
    std::vector<scalar_t> signal_times, signal_values;
    DeltaVector delta_vec =  CoordUtils::MakeDeltaVectorTXYZ(.5, .5, .5, .5);
/*
    std::cout << "Getting current \n";
    SparseCurrentDensity3D current = shower.get_current(delta_vec);
    std::cout << "Starting radio signal calculation \n";
    for(scalar_t cur_t = 1050; cur_t < 1300; cur_t += 1. / sampling_rate) {
        // std::cout << cur_t << "\n";
        scalar_t cur_signal = signal_calc.ComputeSignal(current, cur_t);
        signal_times.push_back(cur_t);
        signal_values.push_back(cur_signal);
    }

    for (int i; i < signal_values.size(); i++) {
        signal_values[i] = signal_values[i] * (1. / constants::epsilon_0 / constants::c / constants::c);
    }
    ExportSignal(signal_times, signal_values, output_path);
    */
    return 0;
}