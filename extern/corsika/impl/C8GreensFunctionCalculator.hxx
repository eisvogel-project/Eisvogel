#include "C8GreensFunctionCalculator.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"

C8CylindricalGreensFunctionCalculator::C8CylindricalGreensFunctionCalculator(float r_max, float z_min, float z_max, float z_ant,
									     std::function<eps_signature> eps,
									     std::function<impulse_response_signature> impulse_response) {

}

void C8CylindricalGreensFunctionCalculator::calculate(std::filesystem::path gf_path,
						      std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir) {

}
