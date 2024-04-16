#include <fstream>
#include <iostream>
#include "Eisvogel/Common.hh"
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

  RZTCoordVector start_coords{0.0f, 0.0f, 0.0f};
  RZTCoordVector end_coords{300.0f, 15.0f, 300.0f};
  RZTVector<std::size_t> num_samples{30u, 2u, 30u};
  
  gf.fill_array<Interpolation::Kernel::Keys>(start_coords, end_coords, num_samples, array);

  std::filesystem::path outpath = "./snapshot.bin";
  std::fstream iofs;  
  iofs.open(outpath, std::ios::out | std::ios::binary);  
  stor::Traits<CylindricalGreensFunction::chunk_t>::serialize_to_numpy(iofs, array);
  iofs.close();
  
  std::cout << "done" << std::endl;
}
