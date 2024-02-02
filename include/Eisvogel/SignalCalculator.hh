#ifndef __SIGNAL_CALCULATOR_HH
#define __SIGNAL_CALCULATOR_HH

#include <filesystem>
#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingField.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"

class SignalCalculator {

public:
  
  SignalCalculator(const std::filesystem::path& geometry_path);
  scalar_t ComputeSignal(const Current0D& track, scalar_t t_sig);
  scalar_t ComputeSignal(const SparseCurrentDensity3D& current_distribution, scalar_t t_sig);
  
private:
  std::filesystem::path m_geometry_path;
  std::shared_ptr<CylindricalWeightingField> m_wf;
};

#endif
