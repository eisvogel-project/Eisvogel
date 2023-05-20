#ifndef __INTEGRATOR_HH
#define __INTEGRATOR_HH

#include <memory>
#include "Common.hh"
#include "WeightingField.hh"
#include "Interpolator.hh"
#include "Current0D.hh"
#include "SparseCurrentDensity3D.hh"

class Integrator {

public:

  void SetGeometry(std::shared_ptr<WeightingField> wf, std::shared_ptr<Kernel> kernel);
  scalar_t integrate(scalar_t t, const Current0D& curr, scalar_t os_factor = 1.0) const;
  scalar_t integrate(scalar_t t, const SparseCurrentDensity3D& current_distribution) const;

private:

  std::shared_ptr<Kernel> m_kernel;
  std::shared_ptr<WeightingField> m_wf;

  using interpolator_t = Interpolator<DenseNDArray, scalar_t, 3>;

  std::unique_ptr<interpolator_t> m_itpl_E_r;
  std::unique_ptr<interpolator_t> m_itpl_E_z;
  std::unique_ptr<interpolator_t> m_itpl_E_phi;
};

#endif
