#ifndef __SIGNAL_CALCULATOR_HH
#define __SIGNAL_CALCULATOR_HH

#include <string>
#include <memory>
#include "Eisvogel/Common.hh"
#include "Eisvogel/Integrator.hh"
#include "Eisvogel/WeightingField.hh"

class SignalCalculator {

public:
   SignalCalculator(const std::string& geometry_path);
   scalar_t ComputeSignal(const Current0D& track, scalar_t t_sig);

private:
   std::string m_geometry_path;
   Integrator m_integrator;
   static std::shared_ptr<WeightingField> load_wf(const std::string& path);
};

#endif
