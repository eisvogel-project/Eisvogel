#pragma once

#include <array>
#include <vector>
#include <filesystem>

class CylindricalGreensFunction;

class C8SignalCalculator {

public:
  C8SignalCalculator(std::filesystem::path gf_path, std::size_t cache_depth = 5);
  
  void calculate(const std::array<float, 4>& track_start_xyzt, const std::array<float, 4>& track_end_xyzt, float track_charge,
		 float t_sig_start, float t_sig_samp, std::size_t num_samples, std::vector<float>& signal);
  
private:
  std::shared_ptr<CylindricalGreensFunction> m_gf;
  
};
