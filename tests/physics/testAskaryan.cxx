#include <filesystem>
#include <cstdlib>

#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"
#include "AnalyticGreensFunctionCalculator.hh"
#include "CSVReader.hh"
#include "testUtils.hh"

void create_greens_function(std::filesystem::path outdir, scalar_t filter_t_peak, unsigned int filter_order) {

  // Build Green's function
  RZCoordVector start_coords{0.0f, -100.0f};
  RZCoordVector end_coords{300.0f, 100.0f};
  scalar_t t_end = 300;
  scalar_t ior = 1.0;
    
  // Sampling parameters
  scalar_t os_factor = 20;
  scalar_t r_min = 0.1;
  
  GreensFunctionCalculator::Analytic::ElectricDipole(outdir, start_coords, end_coords, t_end, ior, filter_t_peak, filter_order, r_min, os_factor);
}

void calculate_signal_eisvogel(std::filesystem::path gf_path, scalar_t b, scalar_t beta,
			       scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples, std::vector<scalar_t>& signal_values) {
  
  // Some parameters are fixed
  scalar_t z0 = 100;
  scalar_t tstart = -250, tend = 250;
  auto charge = [&](scalar_t t){
    return -exp(beta * t / z0);
  };

  // Open Green's function
  SignalCalculator calc(gf_path);

  std::chrono::microseconds total_duration{0};

  for(std::size_t i = 0; i < 2; i++) {
    total_duration = std::chrono::microseconds::zero();
    
    signal_values.resize(num_samples);
    std::fill(signal_values.begin(), signal_values.end(), (scalar_t)0.0f);
    
    // Build trajectory and integrate
    constexpr scalar_t delta_t = 1.0f;
    for(scalar_t cur_t = tstart; cur_t < tend; cur_t += delta_t) {
      
      LineCurrentSegment cur_track(XYZCoordVector{beta * cur_t, 0.0f, b},               // track start position
				   XYZCoordVector{beta * (cur_t + delta_t), 0.0f, b},   // track end position
				   cur_t, cur_t + delta_t, charge(cur_t));
      
      auto start = std::chrono::high_resolution_clock::now();
      calc.AccumulateSignal(cur_track, t_sig_start, t_sig_samp, num_samples, signal_values);
      auto stop = std::chrono::high_resolution_clock::now();
      total_duration += std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    }
  }
  
  std::cout << "Calculation took " << total_duration.count() << "us" << std::endl;
}

void run_test(std::filesystem::path result_path, std::filesystem::path gf_path, scalar_t b, scalar_t beta, float rel_th = 5e-9) {

  // Read in the answer that we are supposed to get
  CSVReader<float> val(result_path);
  std::vector<float> times, signal_result;  
  val.read_column(0, times);
  val.read_column(1, signal_result);
  assert(times.size() == signal_result.size());  
  
  scalar_t t_sig_start = times[0];
  scalar_t t_sig_samp = times[1] - times[0];
  std::size_t num_samples = times.size();
  std::vector<scalar_t> eisvogel_signal_buffer;

  calculate_signal_eisvogel(gf_path, b, beta, t_sig_start, t_sig_samp, num_samples, eisvogel_signal_buffer);
  assert(eisvogel_signal_buffer.size() == signal_result.size());

  std::cout << result_path << std::endl;
  assert(TestUtils::signals_close_match(eisvogel_signal_buffer, signal_result, rel_th, true));
}

int main(void) {

  std::filesystem::path tmpdir = "./testAskaryan/";

  // Prepare Green's function
  scalar_t filter_t_peak = 5.0;
  unsigned int filter_order = 4;  
  create_greens_function(tmpdir, filter_t_peak, filter_order);

  std::filesystem::path validation_dir = std::filesystem::path(std::getenv("PHYSICS_TEST_DATA_DIR")) / "testAskaryan";

  {
    scalar_t b = 10;
    scalar_t beta = 0.9;
    std::filesystem::path result_path = validation_dir / "b_10_beta_0.9_tpeak_5.0_order_4.csv";
    run_test(result_path, tmpdir, b, beta);
  }
  
  // Clean up
  std::filesystem::remove_all(tmpdir);
  
  return 0;
}
