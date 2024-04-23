#include <iostream>
#include <cmath>
#include "Common.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"
#include "Antenna.hh"
#include "Geometry.hh"

#include "CSVReader.hh"

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  if(argc < 3) {
    throw std::runtime_error("Error: need to pass path to output directory and IOR file!");
  }

  std::string gf_path = argv[1];
  std::string ior_path = argv[2];

  // Smooth ice model
  auto eps_smooth = []([[maybe_unused]] scalar_t r, scalar_t z) {
    
    scalar_t z_m = z / 3.0;    
    if(z_m > 0.0) {
      return 1.0;
    }

    scalar_t density = 0.0;

    if(z_m > -14.9) {
      density = 0.917 - 0.594 * std::exp(z_m / 30.8);
    }
    else {
      density = 0.917 - 0.367 * std::exp((z_m + 14.9) / 40.5);
    }

    double eps = std::pow(1 + 0.845 * density, 2.0);
    return eps;
  };

  // Complicated ice model
  CSVReader<float> ior_file(ior_path);
  std::vector<float> depth_m_data;
  std::vector<float> ior_data;
  ior_file.read_column(0, depth_m_data);
  ior_file.read_column(1, ior_data);
  
  auto eps_meas = [&](scalar_t r, scalar_t z) {

    scalar_t z_m = z / 3.0;
    if(z_m > 0.0) {
      return 1.0;      
    }

    if(z_m <= -100.0) {
      return eps_smooth(r, z);
    }
    
    // Now, look up the index of refraction data
    scalar_t depth_m = std::fabs(z_m);

    std::size_t i = 0;
    for(i = 0; i < depth_m_data.size() - 1; i++) {
      if(std::fabs(depth_m_data[i] - depth_m) < 0.1) {
	break;
      }
    }
    
    assert(std::fabs(depth_m_data[i] - depth_m) <= 0.11);
    assert(std::fabs(depth_m_data[i + 1] - depth_m) <= 0.11);
    
    scalar_t ior = ior_data[i] + (ior_data[i+1] - ior_data[i]) / (depth_m_data[i+1] - depth_m_data[i]) * (depth_m - depth_m_data[i]);
    double eps = std::pow(ior, 2.0);
    
    return eps;
  };
  
  auto impulse_response = [](scalar_t t) {
    unsigned int N = 4; // order of filter
    double tp = 2.0; // peaking time of filter
    if(t <= 0) {
      return 0.0;
    }
    return std::pow(t / tp * N, N) * std::exp(-t / tp * N) / (tp * std::exp(std::lgamma(N)));
  };

  // CylinderGeometry geom(20.0, -15.0, 15.0, eps);
  // InfEDipoleAntenna dipole(0.0, 10.0, 0.0, impulse_response);
  // scalar_t t_end = 25.0;

  // CylinderGeometry geom(100.0, -100.0, 100.0, eps);
  // InfEDipoleAntenna dipole(0.0, 10.0, 0.0, impulse_response);
  // scalar_t t_end = 25.0;
  
  CylinderGeometry geom(400, -300, 200, eps_meas);
  InfEDipoleAntenna dipole(0.0, 10.0, -100.0, impulse_response);
  scalar_t t_end = 500;
  
  GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator gfc(geom, dipole, t_end);
  gfc.Calculate(gf_path, "/scratch/midway3/windischhofer/", "/scratch/midway3/windischhofer/");
    
  std::cout << "done" << std::endl;
  
  return 0;
}
