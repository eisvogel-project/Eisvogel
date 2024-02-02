#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Serialization.hh"
#include "IntegratorNew.hh"
#include <iostream>
#include <fstream>

SignalCalculator::SignalCalculator(const std::string& geometry_path) :
  m_geometry_path(geometry_path) {

  // Later, might have to do some dispatch based on the type of the weighting field so that the `SignalCalculator` does not carry a template argument
  m_wf = std::make_shared<CylindricalWeightingField>(m_geometry_path);
}

scalar_t SignalCalculator::ComputeSignal(const Current0D& track, scalar_t t_sig) {
  
  return integrate<CylindricalWeightingField>(*m_wf, t_sig, track);;
}

scalar_t SignalCalculator::ComputeSignal(const SparseCurrentDensity3D& current_distribution, scalar_t t_sig) {

  return integrate<CylindricalWeightingField>(*m_wf, t_sig, current_distribution);  
}
