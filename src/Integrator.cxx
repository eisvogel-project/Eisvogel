#include "Eisvogel/Integrator.hh"
#include "Eisvogel/CoordUtils.hh"

#include <utility>
#include <iostream>

#include "Eisvogel/Trajectory.hh"
#include "Eisvogel/SignalExport.hh"

namespace CU = CoordUtils;

void Integrator::SetGeometry(std::shared_ptr<WeightingField> wf, 
			     std::shared_ptr<Kernel> kernel) {
  m_kernel = kernel;
  m_wf = wf;

  m_itpl_E_r = std::make_unique<interpolator_t>(m_wf -> E_r(), *m_kernel);
  m_itpl_E_z = std::make_unique<interpolator_t>(m_wf -> E_z(), *m_kernel);
  m_itpl_E_phi = std::make_unique<interpolator_t>(m_wf -> E_phi(), *m_kernel);
}

scalar_t Integrator::integrate(scalar_t t, const Current0D& curr, scalar_t os_factor) const {

  // compute velocity vector for each segment
  Trajectory deltas, velocities;
  for(std::size_t pt_ind = 0; pt_ind < curr.number_segments(); pt_ind++) {
    deltas.AddPoint(curr.GetPoint(pt_ind + 1) - curr.GetPoint(pt_ind));
    velocities.AddPoint(deltas(pt_ind) / CU::getT(deltas(pt_ind)));
  }

  DeltaVector wf_sampling_intervals = m_wf -> getSamplingIntervals(); 

  // Main signal integration loop
  scalar_t signal = 0;
  for(std::size_t segment_ind = 0; segment_ind < deltas.size(); segment_ind++) {

    CoordVector& segment_velocity = velocities(segment_ind);
    scalar_t segment_charge = curr.GetCharge(segment_ind);
    
    scalar_t t_step = 1.0 / (1.0 / CU::getT(wf_sampling_intervals) + 
			     std::sqrt(std::pow(CU::getX(segment_velocity), 2) + std::pow(CU::getY(segment_velocity), 2)) / CU::getR(wf_sampling_intervals) + 
			     std::fabs(CU::getZ(segment_velocity)) / CU::getZ(wf_sampling_intervals)
			     );
    t_step /= os_factor;

    scalar_t t_start = CU::getT(curr.GetPoint(segment_ind));
    scalar_t t_end = std::min(t, CU::getT(curr.GetPoint(segment_ind + 1)));

    if(t_end <= t_start) {
      continue;
    }

    // choose final time step so that an integer number of sampling points fits
    const std::size_t number_points = std::ceil((t_end - t_start) / t_step);
    t_step = (t_end - t_start) / number_points;

    // Integrate along segment (bail out early if allowed by causality)
    scalar_t segment_signal = 0;
    scalar_t cur_t = t_start - t_step * m_kernel -> Support();
    for(int step_ind = -m_kernel -> Support(); step_ind <= (int)(number_points + m_kernel -> Support()); step_ind++) {

      CoordVector cur_pos_txyz = curr.GetPoint(segment_ind) + deltas(segment_ind) * (cur_t - t_start) / CU::getT(deltas(segment_ind));
      CoordVector cur_pos_trz = CU::TXYZ_to_TRZ(cur_pos_txyz);
      CoordVector wf_eval_pos = CU::MakeCoordVectorTRZ(t - cur_t, CU::getR(cur_pos_trz), CU::getZ(cur_pos_trz));
      
      CoordVector wf_eval_frac_inds = m_wf -> getFracInds(wf_eval_pos);

      FieldVector wf_rzphi = CU::MakeFieldVectorRZPHI(m_itpl_E_r -> Interpolate(wf_eval_frac_inds),
						      m_itpl_E_z -> Interpolate(wf_eval_frac_inds),
						      m_itpl_E_phi -> Interpolate(wf_eval_frac_inds));

      FieldVector wf_xyz = CU::RZPHI_to_XYZ(wf_rzphi, cur_pos_txyz);

      scalar_t wf_val = CU::getXComponent(wf_xyz) * CU::getX(segment_velocity) +
	CU::getYComponent(wf_xyz) * CU::getY(segment_velocity) +
	CU::getZComponent(wf_xyz) * CU::getZ(segment_velocity);
      
      scalar_t kernel_int = m_kernel -> CDF(number_points - step_ind) - m_kernel -> CDF(-step_ind);
      
      segment_signal += -wf_val * kernel_int;
      cur_t += t_step;
    }
    segment_signal *= t_step * segment_charge;
    signal += segment_signal;
  }

  return signal;
}

scalar_t Integrator::integrate(scalar_t t, const SparseCurrentDensity3D& current_distribution) const {

  scalar_t signal = 0;
  scalar_t volume_element = current_distribution.getVolumeElementTXYZ();

  for(const auto& [cur_pos_txyz, cur_current_density_xyz] : current_distribution) {

    std::cout << "cur_point: t = " << CU::getT(cur_pos_txyz) << ", x = " << CU::getX(cur_pos_txyz) << ", y = " << CU::getY(cur_pos_txyz) << ", z = " << CU::getZ(cur_pos_txyz) << std::endl;
    std::cout << "cur_current: x = " << CU::getXComponent(cur_current_density_xyz) << ", y = " << CU::getYComponent(cur_current_density_xyz) << ", z = " << CU::getZComponent(cur_current_density_xyz) << std::endl;

    scalar_t cur_t = CU::getT(cur_pos_txyz);

    CoordVector cur_pos_trz = CU::TXYZ_to_TRZ(cur_pos_txyz);
    CoordVector wf_eval_pos = CU::MakeCoordVectorTRZ(t - cur_t, CU::getR(cur_pos_trz), CU::getZ(cur_pos_trz));

    CoordVector wf_eval_frac_inds = m_wf -> getFracInds(wf_eval_pos);
    
    FieldVector wf_rzphi = CU::MakeFieldVectorRZPHI(m_itpl_E_r -> Interpolate(wf_eval_frac_inds),
						    m_itpl_E_z -> Interpolate(wf_eval_frac_inds),
						    m_itpl_E_phi -> Interpolate(wf_eval_frac_inds));
    
    FieldVector wf_xyz = CU::RZPHI_to_XYZ(wf_rzphi, cur_pos_txyz);

    scalar_t wf_val = CU::getXComponent(wf_xyz) * CU::getXComponent(cur_current_density_xyz) +
      CU::getYComponent(wf_xyz) * CU::getYComponent(cur_current_density_xyz) +
      CU::getZComponent(wf_xyz) * CU::getZComponent(cur_current_density_xyz);

    signal += -wf_val * volume_element;

  }

  return signal;
}
