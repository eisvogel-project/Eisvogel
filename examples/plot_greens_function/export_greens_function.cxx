#include <fstream>
#include <iostream>
#include "Common.hh"
#include "NDVecArray.hh"
#include "Interpolation.hh"
#include "GreensFunction.hh"

int main(int argc, char* argv[]) {

  if(argc < 3) {
    throw;
  }

  std::filesystem::path gf_path = argv[1];
  std::filesystem::path outpath = argv[2];
  CylindricalGreensFunction gf(gf_path, 5);

  Vector<std::size_t, 3> init_shape(10);
  CylindricalGreensFunction::chunk_t array(init_shape);
  
  RZTVector<std::size_t> num_samples{1000u, 1000u, 1000u};

  RZTCoordVector start_coords = gf.start_coords();
  RZTCoordVector end_coords = gf.end_coords();
  
  std::cout << start_coords << std::endl;
  std::cout << end_coords << std::endl;
  
  // start_coords.r() = 450;
  // start_coords.z() = -400;

  // end_coords.r() = 454;
  // end_coords.z() = -395;

  std::cout << start_coords << std::endl;
  std::cout << end_coords << std::endl;
  
  gf.fill_array<Interpolation::Kernel::Keys>(start_coords, end_coords, num_samples, array);
  
  std::fstream iofs;  
  iofs.open(outpath, std::ios::out | std::ios::binary);  
  stor::Traits<CylindricalGreensFunction::chunk_t>::serialize_to_numpy(iofs, array);
  iofs.close();
  
  std::cout << "done" << std::endl;
}
