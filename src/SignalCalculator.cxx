#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Serialization.hh"
#include "IntegratorNew.hh"
#include <iostream>
#include <fstream>

SignalCalculator::SignalCalculator(const std::string& geometry_path) :
  m_geometry_path(geometry_path), m_integrator()
{ 
  // For now, simply load the weighting field and set it
  // Later, loading the weighting field will be removed here ...
  std::shared_ptr<Kernel> kernel = std::make_shared<KeysCubicInterpolationKernel>();
  std::shared_ptr<DistributedWeightingField> dwf = std::make_shared<DistributedWeightingField>(m_geometry_path);
  m_integrator.SetGeometry(dwf, kernel);

  m_wf = std::make_shared<CylindricalWeightingField>(m_geometry_path);
}

scalar_t SignalCalculator::ComputeSignal(const Current0D& track, scalar_t t_sig) {

  // ... and instead happen here, together with the logic of extracting only a region around
  // the shower (still need to pay attention that we're not generating I/O unecessarily)
  
  // Compute the signal
  // scalar_t sig_old = m_integrator.integrate(t_sig, track);
  scalar_t sig_new = integrate<CylindricalWeightingField>(*m_wf, t_sig, track);

  // std::cout << "old = " << sig_old << "  < --- >   new = " << sig_new << std::endl;
  
  return sig_new;
}

scalar_t SignalCalculator::ComputeSignal(const SparseCurrentDensity3D& current_distribution, scalar_t t_sig) {

  return m_integrator.integrate(t_sig, current_distribution);
}
