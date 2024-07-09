#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include "Common.hh"
#include "CSVReader.hh"
#include "MEEPCylindricalGreensFunctionCalculator.hh"

template <typename T>
std::vector<T> reorder_vector(const std::vector<T>& vec, const std::vector<std::size_t>& inds) {
  std::vector<T> ordered_vec;
  for(std::size_t cur_ind: inds) {
    ordered_vec.push_back(vec[cur_ind]);
  }
  return ordered_vec;
}

// Takes the mapping (argval, funcval) and linearly interpolates the function values to evaluate
// the mapping at `arg_eval`
// Assumes that `argval` is sorted in ascending order
template <typename T>
T lin_interp(const std::vector<T>& argvals, const std::vector<T>& funcvals, const T& arg_eval) {
  assert(argvals.size() == funcvals.size());
  assert(argvals.size() >= 2);

  // Find the two arguments that bracket the requested evaluation point `arg_eval`
  auto find_bracketing_args = [&](const T& arg_a, const T& arg_b) -> bool {
    return (arg_a <= arg_eval) && (arg_b > arg_eval);
  };
  typename std::vector<T>::const_iterator it = std::adjacent_find(argvals.begin(), argvals.end(), find_bracketing_args);

  if(it == argvals.end()) {
    // Out-of-range evaluation requested
    throw std::runtime_error("Error: requested out-of-range evaluation!");
  }

  std::size_t ind_lower = std::distance(argvals.begin(), it);
  std::size_t ind_upper = ind_lower + 1;

  assert(arg_eval >= argvals[ind_lower]);
  assert(arg_eval - argvals[ind_lower] <= argvals[ind_upper] - argvals[ind_lower]);
  
  // Now can do the interpolation
  return funcvals[ind_lower] + (funcvals[ind_upper] - funcvals[ind_lower]) / (argvals[ind_upper] - argvals[ind_lower]) * (arg_eval - argvals[ind_lower]);
}

int main(int argc, char* argv[]) {

  meep::initialize mpi(argc, argv);

  std::filesystem::path gf_path = "/project/avieregg/eisvogel/gf_dipole_summit_pert19_butterworth/";
  std::filesystem::path ior_path = "/home/windischhofer/Eisvogel/applications/dipole/summit_ice_pert19.csv";
  std::filesystem::path impulse_response_path = "/home/windischhofer/Eisvogel/applications/dipole/butterworth_bandpass.csv";
  std::filesystem::path scratch_dir = "/scratch/midway3/windischhofer/";

  scalar_t geom_r_max = 300;  // Radial extent of the simulation domain
  scalar_t geom_z_min = -100;  // Lower z-coordinate of the simulation domain
  scalar_t geom_z_max = 100;  // Upper z-coordinate of the simulation domain
  scalar_t antenna_z = -30;  // z-coordinate of the antenna position (antenna is at r=0 by default)
  scalar_t t_end = 300;  // Last timestep in the simulation

  // Complicated ice model
  CSVReader<float> ior_file(ior_path);
  std::vector<float> z_data;
  std::vector<float> ior_data;
  ior_file.read_column(0, z_data);
  ior_file.read_column(1, ior_data);

  // Ensure increasing z-ordering of the ior data: find the permutation that sorts both arrays
  std::vector<std::size_t> sorter(z_data.size());
  std::iota(sorter.begin(), sorter.end(), 0);
  auto comp = [&](const std::size_t& ind_a, const std::size_t& ind_b) -> bool {
    return z_data[ind_a] < z_data[ind_b];
  };  
  std::sort(sorter.begin(), sorter.end(), comp);

  // Use that permutation to sort them
  z_data = reorder_vector(z_data, sorter);
  ior_data = reorder_vector(ior_data, sorter);

  // Define permittivity based on linearly-interpolated data
  auto eps = [&](scalar_t r, scalar_t z) -> scalar_t {
    (void)r; // The ice model does (for now) not depend on the radial distance from the antenna
        
    scalar_t ior = lin_interp(z_data, ior_data, z);
    double eps = std::pow(ior, 2.0);    
    return eps;
  };

  // Read impulse response from csv file
  CSVReader<float> imp_res_file(impulse_response_path);
  std::vector<float> t_data;
  std::vector<float> imp_res_data;
  imp_res_file.read_column(0, t_data);
  imp_res_file.read_column(1, imp_res_data);

  float imp_res_t_end = t_data.back();
  
  auto impulse_response = [&](scalar_t t) -> scalar_t {

    if((t <= 0) || (t >= imp_res_t_end)) {
      return 0.0;
    }

    scalar_t response = lin_interp(t_data, imp_res_data, t);    
    return response;
  };
  
  CylinderGeometry geom(geom_r_max, geom_z_min, geom_z_max, eps);
  InfEDipoleAntenna dipole(0.0, imp_res_t_end, antenna_z, impulse_response);

  GreensFunctionCalculator::MEEP::CylindricalGreensFunctionCalculator gfc(geom, dipole, t_end);
  gfc.Calculate(gf_path, scratch_dir, scratch_dir);
  
  std::cout << "done" << std::endl;
  
  return 0;
}
