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

  RZTVector<std::size_t> num_samples{2000u, 200u, 2500u};
  gf.fill_array<Interpolation::Kernel::Keys>(gf.start_coords(), gf.end_coords(), num_samples, array);
  
  std::fstream iofs;  
  iofs.open(outpath, std::ios::out | std::ios::binary);  
  stor::Traits<CylindricalGreensFunction::chunk_t>::serialize_to_numpy(iofs, array);
  iofs.close();
  
  std::cout << "done" << std::endl;
}
