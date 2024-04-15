#ifndef __SIGNAL_CALCULATOR_OLD_HH
#define __SIGNAL_CALCULATOR_OLD_HH

#include <filesystem>
#include "Eisvogel/Common.hh"
#include "Eisvogel/WeightingField.hh"
#include "Eisvogel/Current0DOld.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"

class SignalCalculatorOld {

public:
  
  SignalCalculatorOld(const std::filesystem::path& geometry_path);
  scalar_t ComputeSignal(const Current0D& track, scalar_t t_sig);
  scalar_t ComputeSignal(const SparseCurrentDensity3D& current_distribution, scalar_t t_sig);

   /**
     * Calculate the contribution to the antenna signal arising from a particle track.
     *
     * @param track    The particle track.
     * @param ts       The time values at which the signal should be calculated
     * @param signal   Antenna signal buffer in which the contribution is to be accumulated.
     *                 Must exist and have the same length as ts.
     */  
  // void AccumulateSignal(const Current0D& track, std::vector<scalar_t>& ts, std::vector<scalar_t>& signal);
  
private:
  std::filesystem::path m_geometry_path;
  std::shared_ptr<CylindricalWeightingField> m_wf;
};

#endif
