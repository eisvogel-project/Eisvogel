#pragma once

#include <filesystem>
#include <memory>
#include <vector>
#include "Common.hh"

class CylindricalGreensFunction;
class LineCurrentSegment;

class SignalCalculator {

public:

  SignalCalculator(const std::filesystem::path& geometry_path, std::size_t cache_depth = 5);
  void AccumulateSignal(const LineCurrentSegment& seg, scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples, std::vector<scalar_t>& signal);

private:

  std::filesystem::path m_geometry_path;
  std::shared_ptr<CylindricalGreensFunction> m_gf;
  
};
