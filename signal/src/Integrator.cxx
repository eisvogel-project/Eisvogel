#include "Integrator.hh"
#include "CoordUtils.hh"

#include <utility>
#include <iostream>

namespace CU = CoordUtils;

Integrator::Integrator(const WeightingField& wf, const Kernel& kernel) : 
  m_wf(wf), m_itpl_E_r(wf.E_r(), kernel), m_itpl_E_z(wf.E_z(), kernel), m_itpl_E_phi(wf.E_phi(), kernel) { }

scalar_t Integrator::integrate(scalar_t t, const Trajectory& traj, scalar_t os_factor) const {

  // compute velocity vector for each segment
  Trajectory deltas, velocities;
  for(std::size_t pt_ind = 0; pt_ind < traj.size() - 1; pt_ind++) {
    deltas.AddPoint(traj(pt_ind + 1) - traj(pt_ind));
    velocities.AddPoint(deltas(pt_ind) / CU::getT(deltas(pt_ind)));
  }

  // =======================================

  // std::cout << "-------" << std::endl;

  // for(const auto& cur_pt : std::as_const(traj)) {
  //   std::cout << "t = " << CU::getT(cur_pt) << ", x = " << CU::getX(cur_pt) << ", y = " << CU::getY(cur_pt) << ", z = " << CU::getZ(cur_pt) << std::endl;
  // }

  // std::cout << "-------" << std::endl;

  // for(const auto& cur_delta : std::as_const(deltas)) {
  //   std::cout << "delta_t = " << CU::getT(cur_delta) << ", delta_x = " << CU::getX(cur_delta) << ", delta_y = " << CU::getY(cur_delta) << ", delta_z = " << CU::getZ(cur_delta) << std::endl;
  // }  

  // for(const auto& cur_vel : std::as_const(velocities)) {
  //   std::cout << "v_t = " << CU::getT(cur_vel) << ", v_x = " << CU::getX(cur_vel) << ", v_y = " << CU::getY(cur_vel) << ", v_z = " << CU::getZ(cur_vel) << std::endl;
  // }

  // std::cout << "-------" << std::endl;

  // =======================================

  DeltaVector wf_sampling_intervals = m_wf.getSamplingIntervals();
  
  // Main signal integration loop
  scalar_t signal = 0;
  for(std::size_t segment_ind = 0; segment_ind < deltas.size(); segment_ind++) {

    CoordVector& segment_velocity = velocities(segment_ind);
    
    // scalar_t t_step = CU::getT(deltas(segment_ind)) / 500000.0;
    scalar_t t_step = 1.0 / (1.0 / CU::getT(wf_sampling_intervals) + 
			     std::sqrt(std::pow(CU::getX(segment_velocity), 2) + std::pow(CU::getY(segment_velocity), 2)) / CU::getR(wf_sampling_intervals) + 
			     std::fabs(CU::getZ(segment_velocity)) / CU::getZ(wf_sampling_intervals)
			     );
    t_step /= os_factor;

    std::cout << "using t_step = " << t_step << std::endl;

    // 

    scalar_t t_start = CU::getT(traj(segment_ind));
    scalar_t t_end = CU::getT(traj(segment_ind + 1));

    // Integrate along segment (bail out early if allowed by causality)
    for(scalar_t cur_t = t_start; cur_t < std::min(t_end, t); cur_t += t_step) {

      CoordVector cur_pos_txyz = traj(segment_ind) + deltas(segment_ind) * (cur_t - t_start) / CU::getT(deltas(segment_ind));
      CoordVector cur_pos_trz = CU::TXYZ_to_TRZ(cur_pos_txyz);
      CoordVector wf_eval_pos{t - cur_t, CU::getZ(cur_pos_trz), CU::getR(cur_pos_trz)}; // position where to evaluate weighting field (this is like this because of the internal order {t, z, r} in CoordVector ... to be made less error-prone!
      
      CoordVector wf_eval_frac_inds = m_wf.getFracInds(wf_eval_pos);

      // std::cout << "cur_t = " << cur_t << std::endl;
      // std::cout << "t = " << CU::getT(cur_pos_txyz) << ", x = " << CU::getX(cur_pos_txyz) << ", y = " << CU::getY(cur_pos_txyz) << ", z = " << CU::getZ(cur_pos_txyz) << std::endl;
      // std::cout << "t = " << CU::getT(cur_pos_trz) << ", r = " << CU::getR(cur_pos_trz) << ", z = " << CU::getZ(cur_pos_trz) << std::endl;
      // std::cout << "t(eval) = " << CU::getT(wf_eval_pos) << ", r(eval) = " << CU::getR(wf_eval_pos) << ", z(eval) = " << CU::getZ(wf_eval_pos) << std::endl;
      // std::cout << "t_fracind(eval) = " << CU::getT(wf_eval_frac_inds) << ", r_fracind(eval) = " << CU::getR(wf_eval_frac_inds) << ", z_fracind(eval) = " << CU::getZ(wf_eval_frac_inds) << std::endl;

      FieldVector wf_rzphi = CU::MakeFieldVectorRZPHI(m_itpl_E_r.Interpolate(wf_eval_frac_inds),
						      m_itpl_E_z.Interpolate(wf_eval_frac_inds),
						      m_itpl_E_phi.Interpolate(wf_eval_frac_inds));
      
      // std::cout << "E_r = " << CU::getRComponent(wf_rzphi) << ", E_z = " << CU::getZComponent(wf_rzphi) << ", E_phi = " << CU::getPHIComponent(wf_rzphi) << std::endl;

      FieldVector wf_xyz = CU::RZPHI_to_XYZ(wf_rzphi, cur_pos_txyz);

      // std::cout << "E_x = " << CU::getXComponent(wf_xyz) << ", E_y = " << CU::getYComponent(wf_xyz) << ", E_z = " << CU::getZComponent(wf_xyz) << std::endl;

      scalar_t val_step = -t_step * (CU::getXComponent(wf_xyz) * CU::getX(segment_velocity) +
			   CU::getYComponent(wf_xyz) * CU::getY(segment_velocity) +
			   CU::getZComponent(wf_xyz) * CU::getZ(segment_velocity));

      std::cout << "val_step = " << val_step << std::endl;

      signal += val_step;
    }
  }

  return signal;
}
