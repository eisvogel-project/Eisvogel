#include <iostream>
#include <filesystem>
#include "Eisvogel/WeightingFieldCalculator.hh"

WeightingFieldCalculator::WeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna,
						   double courant_factor, double resolution, double pml_width) {
  
  gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  s = std::make_shared<meep::structure>(*gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  f = std::make_shared<meep::fields>(s.get());

  antenna.AddToGeometry(*f, geom);
}

// TODO: later pass buffers to field values that can repeatedly be overwritten as part of chunkloop_data
// (which is a pointer to some custrom struct)
namespace meep {
  void saving_chunkloop(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
			vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
			ivec shift, std::complex<double> shift_phase,
			const symmetry &S, int sn, void *chunkloop_data)
  {   
    // index vectors for start and end of chunk
    ivec isS = S.transform(is, sn) + shift;
    ivec ieS = S.transform(ie, sn) + shift;
    
    // determine rank and shape of this chunk
    int rank = number_of_directions(fc -> gv.dim);
    size_t shape[rank];
    
    int index = 0;
    LOOP_OVER_DIRECTIONS(fc -> gv.dim, d) {
      size_t cur_len = std::max(0, (ie.in_direction(d) - is.in_direction(d)) / 2 + 1);
      shape[index++] = cur_len;
    }
    
    // compute strides
    size_t strides[rank + 1];
    strides[0] = 1;
    for(int i = 0; i < rank; i++) {
      strides[i + 1] = strides[i] * shape[i];
    }
    size_t buflen = strides[rank];
    
    // prepare buffers
    float Ex_buffer[buflen];
    float Ey_buffer[buflen];
    float Ez_buffer[buflen];
    
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
      
      int cur_flat_index = 0;
      index = 0;
      LOOP_OVER_DIRECTIONS(fc -> gv.dim, d) {
	int cur_index = (ichild.in_direction(d) - isS.in_direction(d)) / 2;
	cur_flat_index += cur_index * strides[index];
      }
      
      Ex_buffer[cur_flat_index] = Ex_val.real();
      Ey_buffer[cur_flat_index] = Ey_val.real();
      Ez_buffer[cur_flat_index] = Ez_val.real();        
    } 
    
    // make sure to give enough flexibility to store different fields / field combinations  
    std::string out_dir((char*)chunkloop_data);
    std::string chunk_file("output_chunk_" + std::to_string(ichunk) + ".h5");
    std::string filepath = out_dir + chunk_file;
    
  }
} // end namespace meep
  
void WeightingFieldCalculator::Calculate(double t_end, std::string tmpdir) {

  if(tmpdir.empty()) {
    tmpdir = std::getenv("TMPDIR");
  }
  
  // This is only for cross-checking the geometry for now
  f -> output_hdf5(meep::Dielectric, gv -> surroundings());

  int stepcnt = 0;
  while (f -> time() < t_end) {
    if(meep::am_master()) {
      std::cout << "Simulation time: " << f -> time() << std::endl;
    }
    f -> step();
    std::string data = tmpdir + "/step_" + std::to_string(stepcnt++) + "_";
    // std::string data = "/scratch/midway3/windischhofer/step_" + std::to_string(stepcnt++) + "_";
    f -> loop_in_chunks(meep::saving_chunkloop, (void*)data.c_str(), f -> total_volume());
    // f -> output_hdf5(meep::Ez, gv -> surroundings());
  }
}