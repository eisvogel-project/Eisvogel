#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"
#include "SignalExport.hh"

#include "shower_creator.h"
#include "shower_1D.h"
#include "constants.h"
#include "units.h"

int main(int argc, char* argv[]) {

    scalar_t shower_energy = 1.e18;                  // energy of the shower, in eV
    scalar_t shower_zenith = 90 * units::degree;     // zenith angle of the shower. I recommend leaving this as is
    scalar_t shower_azimuth = 0;                     // azimuth of the shower. I also recommend leaving this
    int is_hadronic = 0;                             // sets the type of the shower. Choose 1 for hardonic or 0 for electromagnetic shower
    int i_shower = 9;                                // index of the shower profile to be chosen from the library. Choose a value between 0 and 9 or set to NULL to have a shower selected at random
    std::array<float, 3> shower_vertex = {-106  , .1, -165};    // position of the shower vertex. If you change this, make sure to also adjust the weighting field accordingly
    scalar_t sampling_rate = 5.;                     // sampling rate that is used when calculating the radio signal

    std::string gf_path;
    std::string output_path;
    if (argc < 2) {
        gf_path = "greens_function";
    } else {
        gf_path = argv[1];
    }
    if (argc < 3) {
        output_path = "shower_radio_emission.csv";
    } else {
        output_path = argv[2];
    }
    SignalCalculator signal_calc(gf_path, 200);
                    
    showers::ShowerCreator shower_creator(std::string(std::getenv("EISVOGELDIR")) + "/extern/shower_profile/shower_file");   
    showers::Shower1D shower = shower_creator.create_shower(shower_vertex, shower_energy, shower_zenith, shower_azimuth, is_hadronic, i_shower);
         
    shower.print_dimensions();

    scalar_t t_sig_start = 1050.0f;
    scalar_t t_sig_end = 1300.0f;
    scalar_t t_sig_samp = 1.0 / sampling_rate;
    std::size_t num_samples = (t_sig_end - t_sig_start) / t_sig_samp;
    std::vector<scalar_t> signal_values(num_samples);

    std::vector<scalar_t> signal_times;
    for(std::size_t sample_ind = 0; sample_ind < num_samples; sample_ind++) {      
      signal_times.push_back(t_sig_start + sample_ind * t_sig_samp);
    }
    assert(signal_times.size() == num_samples);

    std::vector<LineCurrentSegment> tracks;
    shower.fill_tracks(0.2, tracks);

    for(std::size_t i = 0; i < 2; i++) {
    
    std::size_t num_tracks = 0;
    for(auto& cur_track : tracks) {
      std::cout << "have track with start_pos = " << cur_track.start_pos << ", end_pos = " << cur_track.end_pos
		<< ", start_time = " << cur_track.start_time << ", end_time = " << cur_track.end_time << ", charge = " << cur_track.charge << std::endl;

      auto start = std::chrono::high_resolution_clock::now();  
      signal_calc.AccumulateSignal(cur_track, t_sig_start, t_sig_samp, num_samples, signal_values);
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

      num_tracks++;
      
      std::cout << "track took " << duration << std::endl;
    }

    std::cout << "calculated " << num_tracks << " tracks" << std::endl;
    
    }

    ExportSignal(signal_times, signal_values, output_path);
    
    // Current0D current = shower.get_current(0.2);
    // for(scalar_t cur_t = 1050; cur_t < 1300; cur_t += 1. / sampling_rate) {
    //     scalar_t cur_signal = signal_calc.ComputeSignal(current, cur_t);
    //     signal_times.push_back(cur_t);
    //     signal_values.push_back(cur_signal);
    // }
    // // convert to normal units
    // for (int i; i < signal_values.size(); i++) {
    //     signal_values[i] = signal_values[i] * (1. / constants::epsilon_0 / constants::c / constants::c);
    // }
    
    return 0;
}
