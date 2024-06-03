#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"
#include "SignalExport.hh"

int main(int argc, char* argv[]) {

    std::string gf_path;
    std::string output_path;
    if (argc < 2) {
        gf_path = "greens_function";
    } else {
        gf_path = argv[1];
    }

    SignalCalculator signal_calc(gf_path, 200);
                    
    scalar_t t_sig_start = 50.0f;
    scalar_t t_sig_end = 1450.0f;
    scalar_t t_sig_samp = 0.1;
    std::size_t num_samples = (t_sig_end - t_sig_start) / t_sig_samp;
    std::vector<scalar_t> signal_values(num_samples);

    std::vector<scalar_t> signal_times;
    for(std::size_t sample_ind = 0; sample_ind < num_samples; sample_ind++) {      
      signal_times.push_back(t_sig_start + sample_ind * t_sig_samp);
    }
    assert(signal_times.size() == num_samples);

    // fill tracks
    
    std::vector<LineCurrentSegment> tracks;
    scalar_t shower_beta = 1.0;
    scalar_t shower_sigma_x = 10;
    scalar_t shower_t_start = -3 * shower_sigma_x / shower_beta;
    scalar_t shower_t_end = 3 * shower_sigma_x / shower_beta;

    scalar_t shower_z = -240;
    scalar_t shower_y = 0.1;
    scalar_t shower_x_max = 300;
    scalar_t shower_x_start = shower_x_max + 3 * shower_sigma_x;
    scalar_t shower_x_end = shower_x_max - 3 * shower_sigma_x;   

    auto charge = [&](scalar_t cur_t) -> scalar_t {
      scalar_t cur_x = shower_x_start - shower_beta * (cur_t - shower_t_start);
      scalar_t charge = expf(-pow((cur_x - shower_x_max) / shower_sigma_x, 2.0));
      return charge;
    };

    scalar_t delta_t = 0.1;
    for(scalar_t cur_t = shower_t_start; cur_t < shower_t_end; cur_t += delta_t) {

      LineCurrentSegment cur_track(XYZCoordVector{shower_x_start - shower_beta * (cur_t - shower_t_start), shower_y, shower_z},   // track start position
				   XYZCoordVector{shower_x_start - shower_beta * (cur_t + delta_t - shower_t_start), shower_y, shower_z},   // track end position
				   cur_t, cur_t + delta_t, charge(cur_t));
      tracks.push_back(cur_track);
    }

    for(std::size_t i = 0; i < 2; i++) {

      std::size_t total_duration_us = 0;
      
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

	total_duration_us += duration.count();
      }
      
      std::cout << "calculated " << num_tracks << " tracks" << std::endl;
      std::cout << "total shower took " << total_duration_us << "us" << std::endl;
    }

    ExportSignal(signal_times, signal_values, "simple_shower_radio_emission.csv");
    
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
