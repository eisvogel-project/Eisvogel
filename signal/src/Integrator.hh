#ifndef __INTEGRATOR_HH
#define __INTEGRATOR_HH

#include "Common.hh"
#include "WeightingField.hh"
#include "Interpolator.hh"
#include "Trajectory.hh"

class Integrator {

public:

  Integrator(const WeightingField<>& wf, const Kernel& kernel);
  scalar_t integrate(scalar_t t, const Trajectory& traj) const;

private:

  const WeightingField<>& m_wf;

  Interpolator<DenseNDArray, scalar_t, 3> m_itpl_E_r;
  Interpolator<DenseNDArray, scalar_t, 3> m_itpl_E_z;
  Interpolator<DenseNDArray, scalar_t, 3> m_itpl_E_phi;
};

#endif
