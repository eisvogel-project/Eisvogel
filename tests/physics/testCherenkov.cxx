#include <filesystem>
#include <cstdlib>
#include "Common.hh"
#include "SignalCalculator.hh"
#include "AnalyticGreensFunctionCalculator.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"
#include "Antenna.hh"
#include "Geometry.hh"
#include "testUtils.hh"

void create_greens_function_analytic(std::filesystem::path outdir, scalar_t filter_t_peak, unsigned int filter_order) {

  std::filesystem::create_directories(outdir);
  
  // Build Green's function
  RZCoordVector start_coords{0.0f, -10.0f};
  RZCoordVector end_coords{100.0f, 100.0f};
  scalar_t t_end = 100;
  scalar_t ior = 1.0;
    
  // Sampling parameters
  scalar_t os_factor = 50;
  scalar_t r_min = 0.1;
  
  GreensFunctionCalculator::Analytic::ElectricDipole(outdir, start_coords, end_coords, t_end, ior, filter_t_peak, filter_order, r_min, os_factor);
}

void create_greens_function_meep(std::filesystem::path outdir, scalar_t filter_t_peak, unsigned int filter_order,
				 std::filesystem::path scratchdir) {

  auto eps = [](scalar_t r, scalar_t z) {
    (void)r;
    (void)z;

    return 1.0;
  };

  auto impulse_response = [=](scalar_t t) {
    unsigned int N = filter_order; // order of filter
    double tp = filter_t_peak; // peaking time of filter
    if(t <= 0) {
      return 0.0;
    }
    return std::pow(t / tp * N, N) * std::exp(-t / tp * N) / (tp * std::exp(std::lgamma(N)));
  };

  scalar_t t_end = 100;
  CylinderGeometry geom(100, -10, 100, eps);
  InfEDipoleAntenna dipole(0.0, 5 * filter_t_peak, 0.0, impulse_response);

  std::filesystem::create_directories(scratchdir);
  
  GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator gfc(geom, dipole, t_end);
  gfc.Calculate(outdir, scratchdir, scratchdir);
}

void calculate_signal_eisvogel(std::filesystem::path gf_path, scalar_t b, scalar_t beta, scalar_t t_start, scalar_t t_end,
			       scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples, std::vector<scalar_t>& signal_values) {

  scalar_t charge = 1.0;
  
  // Open Green's function
  SignalCalculator calc(gf_path);

  LineCurrentSegment track(XYZCoordVector{beta * t_start, 0.0f, b},    // track start position
			   XYZCoordVector{beta * t_end, 0.0f, b},   // track end position
			   t_start, t_end, charge);

  signal_values.resize(num_samples);
  std::fill(signal_values.begin(), signal_values.end(), (scalar_t)0.0f);
  
  calc.AccumulateSignal(track, t_sig_start, t_sig_samp, num_samples, signal_values);  
}

// Note: very loose tolerance (3e-3) at the moment due to MEEP calculation doing zero-suppression, but analytical calculation not!
void run_test(std::filesystem::path gf_path_analytic, std::filesystem::path gf_path_meep, float rel_th = 3e-3) {

  scalar_t t_start = -5.0;
  scalar_t t_end = 5.0;
  scalar_t b = 70.0;
  scalar_t beta = 1.2;

  scalar_t t_sig_start = 60.0;
  scalar_t t_sig_end = 95.0;
  scalar_t t_sig_samp = 1.0;
  std::size_t num_samples = (t_sig_end - t_sig_start) / t_sig_samp;

  std::vector<scalar_t> gf_analytic_signal_buffer;
  std::vector<scalar_t> gf_meep_signal_buffer;

  calculate_signal_eisvogel(gf_path_analytic, b, beta, t_start, t_end, t_sig_start, t_sig_samp, num_samples, gf_analytic_signal_buffer);
  calculate_signal_eisvogel(gf_path_meep, b, beta, t_start, t_end, t_sig_start, t_sig_samp, num_samples, gf_meep_signal_buffer);

  assert(TestUtils::signals_close_match(gf_analytic_signal_buffer, gf_meep_signal_buffer, rel_th, true));
}

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  scalar_t filter_t_peak = 5.0;
  unsigned int filter_order = 4;

  std::filesystem::path scratch_dir = std::filesystem::temp_directory_path() / "testCherenkov_scratch";    
  std::filesystem::path gf_path_analytic = "./workdir_testCherenkov/gf_analytic";
  std::filesystem::path gf_path_meep = "./workdir_testCherenkov/gf_meep";  
  
  std::cout << "Using scratch_dir = " << scratch_dir << std::endl;
  
  // Prepare Green's functions
  create_greens_function_analytic(gf_path_analytic, filter_t_peak, filter_order);
  create_greens_function_meep(gf_path_meep, filter_t_peak, filter_order, scratch_dir);

  // Run test
  run_test(gf_path_analytic, gf_path_meep);

  std::filesystem::remove_all(scratch_dir);
  std::filesystem::remove_all(gf_path_analytic);
  std::filesystem::remove_all(gf_path_meep);
  
  return 0;
}
