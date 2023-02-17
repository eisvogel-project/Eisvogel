#include "Integrator.hh"

Integrator::Integrator(const WeightingField<>& wf, const Kernel& kernel) : 
  m_wf(wf), m_itpl_E_r(wf.E_r(), kernel), m_itpl_E_z(wf.E_z(), kernel), m_itpl_E_phi(wf.E_phi(), kernel) { }

scalar_t Integrator::integrate(scalar_t t) {
  return 0.0;
}
