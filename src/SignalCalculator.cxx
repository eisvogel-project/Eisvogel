#include "Eisvogel/SignalCalculator.hh"
#include "Eisvogel/Serialization.hh"
#include <iostream>
#include <fstream>

std::shared_ptr<WeightingField> SignalCalculator::load_wf(const std::string& path) {
  std::fstream ifs; 
  ifs.open(path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);
  return std::make_shared<WeightingField>(iser.deserialize<WeightingField>());
}

SignalCalculator::SignalCalculator(const std::string& geometry_path) :
  m_geometry_path(geometry_path), m_integrator()
{ 
  // For now, simply load the weighting field and set it
  // Later, loading the weighting field will be removed here ...
  std::shared_ptr<Kernel> kernel = std::make_shared<KeysCubicInterpolationKernel>();
  m_integrator.SetGeometry(load_wf(m_geometry_path), kernel);
}

scalar_t SignalCalculator::ComputeSignal(const Current0D& track, scalar_t t_sig) {

  // ... and instead happen here, together with the logic of extracting only a region around
  // the shower (still need to pay attention that we're not generating I/O unecessarily)
  
  // Compute the signal
  return m_integrator.integrate(t_sig, track);
}
