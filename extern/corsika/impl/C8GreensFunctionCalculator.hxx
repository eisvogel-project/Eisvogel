#include "C8GreensFunctionCalculator.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"

C8CylindricalGreensFunctionCalculator::C8CylindricalGreensFunctionCalculator(int argc, char* argv[]) {
  m_mpi = std::make_unique<meep::initialize>(argc, argv);
}

C8CylindricalGreensFunctionCalculator::~C8CylindricalGreensFunctionCalculator() = default;

void C8CylindricalGreensFunctionCalculator::calculate(float r_max, float z_min, float z_max, float z_ant, float t_end,
						      std::function<eps_signature> eps,
						      std::function<impulse_response_signature> impulse_response,
						      std::filesystem::path gf_path,
						      std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir) {

  m_geom = std::make_unique<CylinderGeometry>(r_max, z_min, z_max, eps);
  m_ant = std::make_unique<InfEDipoleAntenna>(0.0, t_end, z_ant, impulse_response);
  m_gfc = std::make_unique<GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator>(*m_geom, *m_ant, t_end);
  
  m_gfc -> Calculate(gf_path, local_scratchdir, global_scratchdir);
}
