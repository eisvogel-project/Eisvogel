#include <fstream>
#include <iostream>
#include "Common.hh"
#include "NDVecArray.hh"
#include "Interpolation.hh"
#include "GreensFunction.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw;
  }

  std::string gf_path = argv[1];
  CylindricalGreensFunction gf(gf_path, 5);

  Vector<std::size_t, 3> init_shape(10);
  CylindricalGreensFunction::chunk_t array(init_shape);

  RZTCoordVector start_coords{0.0f, -300.0f, 0.0f};
  RZTCoordVector end_coords{400.0f, 200.0f, 500.0f};
  RZTVector<std::size_t> num_samples{500u, 500u, 300u};

  // RZTCoordVector start_coords{0.0f, -300.0f, 300.0f};
  // RZTCoordVector end_coords{400.0f, 200.0f, 301.0f};
  // RZTVector<std::size_t> num_samples{4000u, 5000u, 2u};
  
  gf.fill_array<Interpolation::Kernel::Keys>(start_coords, end_coords, num_samples, array);

  std::cout << std::format("{}, {}", array[{1u, 1u, 1u}][0], array[{1u, 1u, 1u}][1]) << std::endl;
  
  std::filesystem::path outpath = "./snapshot.bin";
  std::fstream iofs;  
  iofs.open(outpath, std::ios::out | std::ios::binary);  
  stor::Traits<CylindricalGreensFunction::chunk_t>::serialize_to_numpy(iofs, array);
  iofs.close();
  
  std::cout << "done" << std::endl;
}
