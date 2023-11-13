#include <iostream>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  f = std::make_shared<meep::fields>(s.get());

  antenna.AddToGeometry(*f, geom);
}

namespace meep {

void saving_chunkloop(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
		      vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
		      ivec shift, std::complex<double> shift_phase,
		      const symmetry &S, int sn, void *chunkloop_data)
{

  // make sure to give enough flexibility to store different fields / field combinations  
  std::cout << "creating file" << std::endl;
  std::string out_dir((char*)chunkloop_data);  
  std::cout << "out_dir = " << out_dir << std::endl; 

  // index vectors for start and end of chunk
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
  
  double buff[bufsz];

  int metabuff[2];
  metabuff[0] = is.in_direction(R);
  metabuff[1] = is.in_direction(Z);
  
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
    
    // fetch field components at child point
    data.update_values(idx);
    std::complex<double> Ex_val = data.values[0];
    std::complex<double> Ey_val = data.values[1];
    std::complex<double> Ez_val = data.values[2];
    
    ptrdiff_t idx2 = loop_i1 * stride[0] + loop_i2 * stride[2] + loop_i3 * stride[1];
    buff[idx2] = Ez_val.real();
  }

  std::cout << "writing chunk ...";
}

}
  
void WeightingFieldCalculator::Calculate(double t_end) {

  std::string data = "testdata";
  f -> output_hdf5(meep::Dielectric, gv -> surroundings());
    
  while (f -> time() < t_end) {
    if(meep::am_master()) {
      std::cout << "Simulation time: " << f -> time() << std::endl;
    }
    f -> step();
    f -> loop_in_chunks(meep::saving_chunkloop, (void*)data.c_str(), f -> total_volume());
    f -> output_hdf5(meep::Ez, gv -> surroundings());
  }
}
