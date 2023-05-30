#include <meep.hpp>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include "H5Cpp.h"
using namespace meep;

double eps(const vec &p) {
  if (p.z() < 5)
    return 4.0;
  return 1.0;
}

unsigned int fact(unsigned arg) {
  unsigned int retval = 1;

  for(unsigned int cur = 1; cur <= arg; cur++) {
    retval *= cur;
  }
  
  return retval;
}

std::complex<double> srcfunc(double t, void*) {

  unsigned int order = 4;
  double tp = 2.0;
  
  if(t > 0) {
    double retval = 1.0 / (tp * fact(order - 1)) * std::pow(t * (double)order / tp, (double)order) * std::exp(-t * (double)order / tp);
    return retval;
  }
  return 0.0;
}

void my_chunkloop(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
                  vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
                  ivec shift, std::complex<double> shift_phase,
		  const symmetry &S, int sn, void *chunkloop_data)
{

  std::cout << "creating file" << std::endl;

  H5::H5File file("output_chunk_" + std::to_string(ichunk) + ".h5", H5F_ACC_TRUNC);

  std::ofstream asciifile;
  asciifile.open("output_chunk_" + std::to_string(ichunk) + ".txt");
  
  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;
  
  size_t bufsz = 1;
  size_t dims[2];
  size_t stride[3] = {0, 1, 0};
  int rank = 0;
  LOOP_OVER_DIRECTIONS(fc->gv.dim, d) {
    size_t cur_len = std::max(0, (ie.in_direction(d) - is.in_direction(d)) / 2 + 1);
    bufsz *= cur_len;
    dims[rank++] = cur_len;
  }
  stride[2] = dims[0];

  std::cout << "bufsz = " << bufsz << std::endl;  
  std::cout << "rank = " << rank << std::endl;
  std::cout << "dims[0] = " << dims[0] << std::endl;
  std::cout << "dims[1] = " << dims[1] << std::endl;

  std::cout << "stride[0] = " << stride[0] << std::endl;
  std::cout << "stride[1] = " << stride[1] << std::endl;
  std::cout << "stride[2] = " << stride[2] << std::endl;

  hsize_t dimsf[2];
  dimsf[0] = dims[0];
  dimsf[1] = dims[1];
  H5::DataSpace dataspace(2, dimsf);

  hsize_t metadims[1] = {2};
  H5::DataSpace metadataspace(1, metadims);
  
  H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
  H5::IntType metadatatype(H5::PredType::NATIVE_INT);
  
  double buff[bufsz];

  int metabuff[2];
  metabuff[0] = is.in_direction(R);
  metabuff[1] = is.in_direction(Z);

  H5::DataSet metadata = file.createDataSet("chunk_start", metadatatype, metadataspace);
  H5::DataSet dataset = file.createDataSet("Ez", datatype, dataspace);  
  
  // some preliminary setup
  vec rshift(shift * (0.5*fc->gv.inva));  // shift into unit cell for PBC geometries
  
  // prepare the list of field components to fetch at each grid point
  component components[] = {Ex, Ey, Ez};
  chunkloop_field_components data(fc, cgrid, shift_phase, S, sn, 3, components);
  
  // loop over all grid points in chunk
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    // get grid indices and coordinates of parent point
    IVEC_LOOP_ILOC(fc->gv, iparent);  // grid indices
    IVEC_LOOP_LOC(fc->gv, rparent);   // cartesian coordinates
    
    // apply symmetry transform to get grid indices and coordinates of child point
    ivec ichild = S.transform(iparent, sn) + shift;
    vec rchild = S.transform(rparent, sn) + rshift;

    std::cout << "chunk_number = " << ichunk << std::endl;
    std::cout << "r_ind = " << ichild.r() << " z_ind = " << ichild.z() << std::endl;
    std::cout << "r = " << rchild.r() << " z = " << rchild.z() << std::endl;
    
    // fetch field components at child point
    data.update_values(idx);
    std::complex<double> Ex_val = data.values[0];
    std::complex<double> Ey_val = data.values[1];
    std::complex<double> Ez_val = data.values[2];

    std::cout << "loop_i1 = " << loop_i1 << std::endl;
    std::cout << "loop_i2 = " << loop_i2 << std::endl;
    std::cout << "loop_i3 = " << loop_i3 << std::endl;
    
    ptrdiff_t idx2 = loop_i1 * stride[0] + loop_i2 * stride[2] + loop_i3 * stride[1];
    std::cout << "buff_idx (bufsz) = " << idx2 << " (" << bufsz << ")" << std::endl;
    std::cout << Ez_val.real() << std::endl;
    buff[idx2] = Ez_val.real();

    std::cout << "- - -" << std::endl;

    asciifile << Ez_val.real() << " " << std::endl;
  }

  std::cout << "writing chunk ...";

  dataset.write(buff, H5::PredType::NATIVE_DOUBLE);
  metadata.write(metabuff, H5::PredType::NATIVE_INT);
  file.close();
  
  asciifile.close();
}

int main(int argc, char **argv) {
  
  initialize mpi(argc, argv);
  double resolution = 2;
  grid_volume gv = volcyl(10, 10, resolution); // (r, z, resolution)
  
  structure s(gv, eps, pml(1.0), identity(), 0, 0.5);
  fields f(&s);

  f.output_hdf5(Dielectric, gv.surroundings());

  custom_src_time src(srcfunc, NULL);
  f.add_point_source(Ez, src, veccyl(0.0, 3.0));

  //for(double cur_t = 2.0; cur_t <= 40; cur_t += 2) {
  double cur_t = 10;
  while (f.time() < cur_t) {
    f.step();
  }
    
  f.loop_in_chunks(my_chunkloop, NULL, f.total_volume());
  
  f.output_hdf5(Ez, gv.surroundings());
    //}
  
  return 0;
}
