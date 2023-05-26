#include <meep.hpp>
#include <cmath>
#include <iostream>
using namespace meep;

double eps(const vec &p) {
  if (p.z() < 10)
    return 2.0;
  return 1.0;
}

std::complex<double> srcfunc(double t, void*) {
  if(t > 20) {
    return 0.0;
  }
  else {
    return std::exp(-50 * std::pow(t - 10, 2));
  }
  return 0.0;
}

void my_chunkloop(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
                  vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
                  ivec shift, std::complex<double> shift_phase,
		  const symmetry &S, int sn, void *chunkloop_data)
{

  std::cout << "creating file" << std::endl;
  h5file* outfile = new h5file(("output_chunk_" + std::to_string(ichunk) + ".h5").c_str(), h5file::WRITE, false, true);
  
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
  
  outfile->create_or_extend_data("Ez", rank, dims, false, false);

  double buff[bufsz];
  
  // some preliminary setup
  vec rshift(shift * (0.5*fc->gv.inva));  // shift into unit cell for PBC geometries
  
  // prepare the list of field components to fetch at each grid point
  component components[] = {Ex, Ey, Ez};
  chunkloop_field_components data(fc, cgrid, shift_phase, S, sn, 2, components);
  
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
    std::cout << Ez_val << std::endl;

    std::cout << "loop_i1 = " << loop_i1 << std::endl;
    std::cout << "loop_i2 = " << loop_i2 << std::endl;
    std::cout << "loop_i3 = " << loop_i3 << std::endl;
    
    ptrdiff_t idx2 = loop_i1 * stride[0] + loop_i2 * stride[2] + loop_i3 * stride[1];
    std::cout << "buff_idx (bufsz) = " << idx2 << " (" << bufsz << ")" << std::endl;
    buff[idx2] = Ez_val.real();

    std::cout << "- - -" << std::endl;
  }

  std::cout << "writing chunk ...";
  outfile->write("Ez", rank, dims, buff, false);
  std::cout << " done!" << std::endl;
  
  delete outfile;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  double resolution = 3;
  grid_volume gv = volcyl(10, 10, resolution); // (r, z, resolution)
  
  structure s(gv, eps, pml(1.0), identity(), 0, 0.5);
  fields f(&s);

  f.output_hdf5(Dielectric, gv.surroundings());

  custom_src_time src(srcfunc, NULL);
  f.add_point_source(Ez, src, veccyl(0.0, 8.0));
  
  while (f.time() < 20) {
     f.step();
  }

  f.loop_in_chunks(my_chunkloop, NULL, f.total_volume());
  
  f.output_hdf5(Ez, gv.surroundings());

  return 0;
}
