#ifndef __INTEGRATOR_NEW_HH
#define __INTEGRATOR_NEW_HH

#include <memory>
#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingField.hh"
#include "Eisvogel/Current0D.hh"
#include "Eisvogel/Kernels.hh"

template <class WeightingFieldT, typename KernelT = KeysCubicInterpolationKernelNew>
scalar_t integrate(WeightingFieldT& wf, scalar_t t, const Current0D& curr, scalar_t os_factor = 1.0);

template <class WeightingFieldT, typename KernelT = KeysCubicInterpolationKernelNew>
scalar_t integrate(WeightingFieldT& wf, scalar_t t, const SparseCurrentDensity3D& current_distribution);

#include "IntegratorNew.hxx"

#endif
