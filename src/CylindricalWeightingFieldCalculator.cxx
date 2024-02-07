#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>
#include "Eisvogel/CylindricalWeightingFieldCalculator.hh"
#include "Eisvogel/WeightingField.hh"
#include "FieldStorage.hh"
#include "Eisvogel/CoordUtils.hh"

// For now, this only handles geometries with cylindrical symmetry
CylindricalWeightingFieldCalculator::CylindricalWeightingFieldCalculator(CylinderGeometry& geom, const Antenna& antenna, scalar_t t_end,
									 double courant_factor, double resolution, double pml_width) {

  m_t_end = t_end;
  m_start_coords = std::make_shared<CoordVector>(CoordUtils::MakeCoordVectorTRZ(0.0, 0.0, geom.GetZMin()));
  m_end_coords = std::make_shared<CoordVector>(CoordUtils::MakeCoordVectorTRZ(t_end, geom.GetRMax(), geom.GetZMax()));
  
  m_gv = std::make_shared<meep::grid_volume>(meep::volcyl(geom.GetRMax(), geom.GetZMax() - geom.GetZMin(), resolution));
  m_s = std::make_shared<meep::structure>(*m_gv, geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
  m_f = std::make_shared<meep::fields>(m_s.get());

  antenna.AddToGeometry(*m_f, geom);
}

struct ChunkloopData {

  ChunkloopData(std::size_t ind_t, RZFieldStorage& fstor) :
    ind_t(ind_t), fstor(fstor) { }
  
  std::size_t ind_t;
  RZFieldStorage& fstor;
  
};

namespace meep {
  void saving_chunkloop(fields_chunk* fc, int ichunk, component cgrid, ivec is, ivec ie,
			vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
			ivec shift, std::complex<double> shift_phase,
			const symmetry& S, int sn, void* cld)
  {

    ChunkloopData* chunkloop_data = static_cast<ChunkloopData*>(cld);
    
    // index vectors for start and end of chunk
    ivec isS = S.transform(is, sn) + shift;
    ivec ieS = S.transform(ie, sn) + shift;
    
    // determine rank and shape of this chunk
    std::size_t rank = number_of_directions(fc -> gv.dim);
    constexpr std::size_t required_rank = 2;
										 
    if(rank != required_rank) {
      throw std::runtime_error("Error: Expected 2d situation!");
    }
										 
    std::size_t shape[required_rank] = {0};
    std::size_t index = 0;
    LOOP_OVER_DIRECTIONS(fc -> gv.dim, d) {
      std::size_t cur_len = std::max(0, (ie.in_direction(d) - is.in_direction(d)) / 2 + 1);
      shape[index++] = cur_len;
    }

    // Prepare chunk buffers
    std::size_t pts_t = 1, pts_r = shape[1], pts_z = shape[0];
    ScalarField3D<scalar_t> chunk_buffer_E_r({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> chunk_buffer_E_z({pts_t, pts_z, pts_r}, 0.0);

    assert(is.z() >= 0);
    assert(is.r() >= 0);
    
    IndexVector chunk_start_inds = {
      chunkloop_data -> ind_t,
      (std::size_t)((is.z() - 1) / 2),
      (std::size_t)((is.r() - 1) / 2)
    };
    
    // std::cout << "shape = " << process_rank << " --> shape: r = " << pts_r << ", z = " << pts_z << std::endl;
        
    // some preliminary setup
    vec rshift(shift * (0.5*fc->gv.inva));  // shift into unit cell for PBC geometries
    
    // prepare the list of field components to fetch at each grid point
    component components[] = {Ez, Er};
    chunkloop_field_components data(fc, cgrid, shift_phase, S, sn, 2, components);
    
    // loop over all grid points in chunk
    // TODO: in principle, need only iterate over all parent points and store those, no need to separately store all points on a symmetry orbit
    LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
      
      // get grid indices and coordinates of parent point
      IVEC_LOOP_ILOC(fc->gv, iparent);  // grid indices
      IVEC_LOOP_LOC(fc->gv, rparent);   // cartesian coordinates
      
      // apply symmetry transform to get grid indices and coordinates of child point
      ivec ichild = S.transform(iparent, sn) + shift;
      vec rchild = S.transform(rparent, sn) + rshift;	    

      IndexVector global_ind = {
	chunkloop_data -> ind_t,
	(std::size_t)((ichild.z() - 1) / 2),
	(std::size_t)((ichild.r() - 1) / 2)
      };
      
      // fetch field components at child point
      data.update_values(idx);
      double E_z_val = data.values[0].real();
      double E_r_val = data.values[1].real();
      
      IndexVector chunk_ind = global_ind - chunk_start_inds;
      
      chunk_buffer_E_r(chunk_ind) = E_r_val;
      chunk_buffer_E_z(chunk_ind) = E_z_val;
    }

    auto to_keep = [](float value) -> bool {
      return value < 1e-6;
    };
    
    SparseScalarField3D<scalar_t> sparse_chunk_buffer_E_r = SparseScalarField3D<scalar_t>::FromDense(chunk_buffer_E_r, to_keep, 0.0);
    SparseScalarField3D<scalar_t> sparse_chunk_buffer_E_z = SparseScalarField3D<scalar_t>::FromDense(chunk_buffer_E_z, to_keep, 0.0);
    
    chunkloop_data -> fstor.RegisterChunk(sparse_chunk_buffer_E_r, sparse_chunk_buffer_E_z, chunk_start_inds);       
  }
  
} // end namespace meep
  
void CylindricalWeightingFieldCalculator::Calculate(std::filesystem::path outdir, std::filesystem::path tmpdir) {

  // TODO: this will get a lot easier once we can have all three components together in the same array
  std::filesystem::path outdir_Er = outdir / "E_r";
  std::filesystem::path outdir_Ez = outdir / "E_z";
  
  if(meep::am_master()) {
    
    // Prepare output directory
    if(!std::filesystem::exists(outdir)) {
      std::filesystem::create_directory(outdir);
    }

    std::filesystem::create_directory(outdir_Er);
    std::filesystem::create_directory(outdir_Ez);
  }

  meep::all_wait();
  
  // Make sure to have a valid temporary directory to store the results as we go along
  if(tmpdir.empty()) {
    char tmpdir_template[] = "/tmp/eisvogel.XXXXXX";
    tmpdir = std::string(mkdtemp(tmpdir_template));
  }

  std::cout << "Using tmpdir = " << tmpdir << std::endl;
  
  // This is only for cross-checking the geometry for now
  // f -> output_hdf5(meep::Dielectric, gv -> surroundings());

  std::shared_ptr<RZFieldStorage> fstor = std::make_shared<RZFieldStorage>(tmpdir, 10);
  ChunkloopData cld(0, *fstor);

  std::size_t stepcnt = 0;
  for(double cur_t = 0.0; cur_t <= m_t_end; cur_t += 0.1) {

    // Time-step the fields
    while (m_f -> time() < cur_t) {
      m_f -> step();
    }
    
    if(meep::am_master()) {
      std::cout << "Simulation time: " << m_f -> time() << std::endl;
    }

    cld.ind_t = stepcnt++;
    m_f -> loop_in_chunks(meep::saving_chunkloop, static_cast<void*>(&cld), m_f -> total_volume());
  }
  
  // TODO: again, will get better once the three separate arrays are gone
  std::cout << "moving output into final location ...";
  std::filesystem::copy(tmpdir / "E_r", outdir_Er, std::filesystem::copy_options::recursive);
  std::filesystem::copy(tmpdir / "E_z", outdir_Ez, std::filesystem::copy_options::recursive);
  std::cout << " done!" << std::endl;

  meep::all_wait();
  
  if(meep::am_master()) {
    std::shared_ptr<CylindricalWeightingField> cwf = std::make_shared<CylindricalWeightingField>(outdir, *m_start_coords, *m_end_coords);
    cwf -> MakeMetadataPersistent();
  }
}
