#include "Eisvogel/SignalCalculatorOld.hh"
#include "IntegratorOld.hh"

SignalCalculatorOld::SignalCalculatorOld(const std::filesystem::path& geometry_path) :
  m_geometry_path(geometry_path) {

  // Later, might have to do some dispatch based on the type of the weighting field so that the `SignalCalculatorOld` does not carry a template argument
  m_wf = std::make_shared<CylindricalWeightingField>(m_geometry_path);
}

scalar_t SignalCalculatorOld::ComputeSignal(const Current0D& track, scalar_t t_sig) {  
  return integrate<CylindricalWeightingField>(*m_wf, t_sig, track);
}

scalar_t SignalCalculatorOld::ComputeSignal(const SparseCurrentDensity3D& current_distribution, scalar_t t_sig) {
  return integrate<CylindricalWeightingField>(*m_wf, t_sig, current_distribution);  
}

// void SignalCalculatorOld::AccumulateSignal(const Current0D& track, std::vector<scalar_t>& ts, std::vector<scalar_t>& signal) {
//   return accumulate_integral<CylindricalWeightingField>(*m_wf, ts, track, signal);
// }
