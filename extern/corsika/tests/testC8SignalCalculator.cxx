#include <iostream>
#include <filesystem>
#include "C8SignalCalculator.hh"

int main(void) {

  std::filesystem::path gf_path = "/home/windischhofer/data/windischhofer/eisvogel/gf_dipole_smooth_ice_deep_z_100_r_1000";
  std::shared_ptr<C8SignalCalculator> calc = std::make_shared<C8SignalCalculator>(gf_path);

  using scalar_t = float;
  
  // --------------
  // debug values
  std::size_t number_samples = 100;
  std::array<scalar_t, 4> track_start_xyzt{
    -20.0f, 0.0f, -10.0f, 10.0f
  };
  std::array<float, 4> track_end_xyzt{
    -10.0f, 0.0f, -10.0f, 30.0f
  };
  scalar_t track_charge = 1.0f;
  scalar_t t_sig_start_ns = 10.0f;
  scalar_t t_sig_samp_ns = 1.0f;
  // --------------
  
  // Calculate the signal
  std::vector<scalar_t> sigbuf(number_samples, 0.0f);
  calc -> calculate(track_start_xyzt, track_end_xyzt, track_charge, t_sig_start_ns, t_sig_samp_ns,
		      number_samples, sigbuf);

}
