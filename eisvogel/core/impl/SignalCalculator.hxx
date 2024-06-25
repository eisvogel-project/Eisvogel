#include "Current.hh"
#include "Common.hh"
#include "Interpolation.hh"
#include "GreensFunction.hh"

SignalCalculator::SignalCalculator(const std::filesystem::path& geometry_path,
                                   std::size_t cache_depth)
    : m_geometry_path(geometry_path)
{

  // Later, need to do some dispatch based on the type of the weighting field so
  // that the `SignalCalculator` does not need to carry a template argument
  m_gf =
      std::make_shared<CylindricalGreensFunction>(m_geometry_path, cache_depth);
}

void SignalCalculator::AccumulateSignal(const LineCurrentSegment& seg,
                                        scalar_t t_sig_start,
                                        scalar_t t_sig_samp,
                                        std::size_t num_samples,
                                        std::vector<scalar_t>& signal)
{

  m_gf->apply_accumulate<Interpolation::Kernel::Keys>(
      seg, t_sig_start, t_sig_samp, num_samples, signal);
}
