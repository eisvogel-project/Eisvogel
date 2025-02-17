#pragma once
#define NDEBUG

#include <array>
#include <vector>
#include <filesystem>

class CylindricalGreensFunction;

class C8SignalCalculator {
  
public:
  C8SignalCalculator(std::filesystem::path gf_path, std::size_t cache_depth = 5);
  ~C8SignalCalculator();

  void get_antenna_location(std::array<float, 3>& ant_xyz);
    
  void calculate(const std::array<float, 4>& track_start_xyzt, const std::array<float, 4>& track_end_xyzt, float track_charge,
		 float t_sig_start, float t_sig_samp, std::size_t num_samples, std::vector<double>& signal, float weight = 1.0);
  
private:

  std::unique_ptr<CylindricalGreensFunction> m_gf;
  
};
