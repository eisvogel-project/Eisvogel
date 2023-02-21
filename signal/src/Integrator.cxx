#include "Integrator.hh"

#include <utility>
#include <iostream>

Integrator::Integrator(const WeightingField<>& wf, const Kernel& kernel) : 
  m_wf(wf), m_itpl_E_r(wf.E_r(), kernel), m_itpl_E_z(wf.E_z(), kernel), m_itpl_E_phi(wf.E_phi(), kernel) { }

scalar_t Integrator::integrate(scalar_t t, const Trajectory& traj) const {

  for(const auto& cur_pt : std::as_const(traj)) {
    for(const auto& coord : std::as_const(cur_pt)) {
      std::cout << coord << " ";
    }
    std::cout << std::endl;
  }

  return 0.0;
}
