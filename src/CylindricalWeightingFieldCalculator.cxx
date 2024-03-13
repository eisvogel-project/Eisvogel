#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Eisvogel/CylindricalWeightingFieldCalculator.hh"
#include "Eisvogel/WeightingField.hh"
#include "FieldStorage.hh"
#include "Eisvogel/IteratorUtils.hh"
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

template <typename KeyT, std::size_t dims>
class FieldStatisticsTracker {

public:

  using shape_t = IndexVector;
  using region_t = DenseNDArray<scalar_t, dims>;
  
  FieldStatisticsTracker(std::size_t subsampling) :
    m_subsampling(subsampling) { }

  void AddRegion(KeyT& region_key, const shape_t& region_shape, scalar_t init_value = 0.0) {
    shape_t subsampled_shape = ToSubsampledShape(region_shape);

    // --- this is just a crutch for now until we have fixed-size vectors
    std::array<std::size_t, dims> subsampled_shape_crutch;
    std::copy(std::begin(subsampled_shape), std::end(subsampled_shape), std::begin(subsampled_shape_crutch));
    // ---------
    
    m_max_data.emplace(std::make_pair(region_key, region_t(subsampled_shape_crutch, init_value)));

    std::cout << "for chunk with key = " << region_key << std::endl;
    std::cout << "constructed max tracker with shape" << std::endl;
    subsampled_shape.print();
    std::cout << "for original chunk with shape" << std::endl;
    region_shape.print();
  }

  void UpdateStatisticsForRegion(KeyT& region_key, const region_t& region_data) {

    auto region_el = m_max_data.find(region_key);
    if(region_el == m_max_data.end()) {
      throw std::runtime_error("Error: requested region not found!");
    }
    region_t& region_max_data = region_el -> second;
    
    IndexVector start_inds(dims, 0);
    IndexVector end_inds = region_data.shape();     
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
      
      IndexVector cur_ind = cnt.index();
      IndexVector cur_subsampled_ind = ToSubsampledInd(cur_ind);            
      scalar_t cur_val = region_data(cur_ind);

      // have a new local field maximum	
      if(cur_val > region_max_data(cur_subsampled_ind)) {       
	region_max_data(cur_subsampled_ind) = cur_val;
      }
    }    
  }

  scalar_t GetMaxLocal(KeyT& region_key, IndexVector& region_inds) {
    auto region_el = m_max_data.find(region_key);
    if(region_el == m_max_data.end()) {
      throw std::runtime_error("Error: requested region not found!");
    }
    IndexVector cur_subsampled_ind = ToSubsampledInd(region_inds);
    return (region_el -> second)(cur_subsampled_ind);
  }  

private:

  IndexVector ToSubsampledShape(const IndexVector& shape) {
    // divide + round up with integer divisions
    return (shape + m_subsampling - 1) / m_subsampling;
  }
  
  IndexVector ToSubsampledInd(const IndexVector& ind) {
    return ind / m_subsampling;
  }
  
private:

  std::size_t m_subsampling;  
  std::map<KeyT, region_t> m_max_data;
  
};

template <std::size_t dims>
struct SimulationChunkMetadata {

  using chunk_shape_t = typename DenseNDArray<scalar_t, dims>::shape_t;
  chunk_shape_t chunk_shape;
  
};

// t = const chunks are 2D for a cylindrical geometry, and Meep indexes its chunks
// through an `int`
using FieldChunkStatisticsTracker2D = FieldStatisticsTracker<int, 2>;
using MeepChunkMetadata = SimulationChunkMetadata<2>;

struct ChunkloopData {

  ChunkloopData(std::size_t ind_t, RZFieldStorage& fstor, FieldChunkStatisticsTracker2D& fstats, scalar_t dynamic_range, scalar_t abs_min_field) :
    ind_t(ind_t), fstor(fstor), fstats(fstats), dynamic_range(dynamic_range), abs_min_field(abs_min_field) { }

  std::size_t ind_t;
  RZFieldStorage& fstor;
  FieldChunkStatisticsTracker2D& fstats;
  scalar_t dynamic_range;
  scalar_t abs_min_field;
  std::map<int, MeepChunkMetadata> chunk_meta;  
  
};

namespace meep {

  void eisvogel_setup_chunkloop(fields_chunk* fc, int ichunk, component cgrid, ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
				ivec shift, std::complex<double> shift_phase,
				const symmetry& S, int sn, void* cld) {
    
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

    // Record chunk metadata
    MeepChunkMetadata this_chunk_meta;
    this_chunk_meta.chunk_shape = {shape[0], shape[1]};        
    chunkloop_data -> chunk_meta[ichunk] = this_chunk_meta;

    // Setup field statistics tracker for this chunk
    chunkloop_data -> fstats.AddRegion(ichunk, this_chunk_meta.chunk_shape);
  }
  
  void eisvogel_saving_chunkloop(fields_chunk* fc, int ichunk, component cgrid, ivec is, ivec ie,
				 vec s0, vec s1, vec e0, vec e1, double dV0, double dV1,
				 ivec shift, std::complex<double> shift_phase,
				 const symmetry& S, int sn, void* cld) {
    
    ChunkloopData* chunkloop_data = static_cast<ChunkloopData*>(cld);
    
    // // index vectors for start and end of chunk
    // ivec isS = S.transform(is, sn) + shift;
    // ivec ieS = S.transform(ie, sn) + shift;
    
    // // determine rank and shape of this chunk
    // std::size_t rank = number_of_directions(fc -> gv.dim);
    // constexpr std::size_t required_rank = 2;
										 
    // if(rank != required_rank) {
    //   throw std::runtime_error("Error: Expected 2d situation!");
    // }
										 
    // std::size_t shape[required_rank] = {0};
    // std::size_t index = 0;
    // LOOP_OVER_DIRECTIONS(fc -> gv.dim, d) {
    //   std::size_t cur_len = std::max(0, (ie.in_direction(d) - is.in_direction(d)) / 2 + 1);
    //   shape[index++] = cur_len;
    // }

    // ------------------
    // Prepare chunk buffers (different chunks will have different sizes; need to allocate a buffer for each chunk)
    // TODO: can have a "scanning" chunkloop that is run once at the very beginning to figure out how many chunks there
    // are and what their dimensions are, then allocate static buffers, and then use those for saving
    // (i.e. just like done for the FieldStatisticsTracker at the moment)
    // ------------------
    MeepChunkMetadata::chunk_shape_t chunk_shape = chunkloop_data -> chunk_meta[ichunk].chunk_shape;
    // std::size_t pts_t = 1, pts_r = shape[1], pts_z = shape[0];
    std::size_t pts_t = 1, pts_r = chunk_shape[1], pts_z = chunk_shape[0];
    ScalarField3D<scalar_t> chunk_buffer_E_r({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> chunk_buffer_E_z({pts_t, pts_z, pts_r}, 0.0);
    
    ScalarField2D<scalar_t> chunk_buffer_Evecnorm_vals({pts_z, pts_r}, 0.0);

    // std::cout << "simulation: have chunk " << ichunk << " with chunk_shape[0] = " << chunk_shape[0] << ", chunk_shape[1] = " << chunk_shape[1] << std::endl;

    assert(is.z() >= 0);
    assert(is.r() >= 0);
    
    IndexVector chunk_start_inds = {
      chunkloop_data -> ind_t,
      (std::size_t)((is.z() - 1) / 2),
      (std::size_t)((is.r() - 1) / 2)
    };
    IndexVector chunk_start_inds_tslice = {
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
      IndexVector global_ind_tslice = {
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

      // use field vector norm to keep track of field maximum
      IndexVector chunk_ind_tslice = global_ind_tslice - chunk_start_inds_tslice;
      chunk_buffer_Evecnorm_vals(chunk_ind_tslice) = std::sqrt(std::pow(E_r_val, 2) + std::pow(E_z_val, 2));
    }

    chunkloop_data -> fstats.UpdateStatisticsForRegion(ichunk, chunk_buffer_Evecnorm_vals);

    // keep track of smallest surviving nonzero absolute values
    scalar_t min_abs_E_r = std::numeric_limits<scalar_t>::max();
    scalar_t min_abs_E_z = std::numeric_limits<scalar_t>::max();
    
    // Restrict dynamic range of output fields based on field statistics
    std::size_t num_truncations = 0;
    IndexVector start_inds(3, 0);
    IndexVector end_inds = chunk_buffer_E_r.shape();     
    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {
      IndexVector cur_ind = cnt.index();
      IndexVector cur_slice_ind = {
	CoordUtils::getZInd(cur_ind),
	CoordUtils::getRInd(cur_ind)
      };
      scalar_t max_field = chunkloop_data -> fstats.GetMaxLocal(ichunk, cur_slice_ind);
	
      if(max_field > chunkloop_data -> abs_min_field) {
	
	scalar_t cur_E_r_abs = std::fabs(chunk_buffer_E_r(cur_ind));
	scalar_t cur_E_z_abs = std::fabs(chunk_buffer_E_z(cur_ind));
	
	if(cur_E_r_abs < max_field / chunkloop_data -> dynamic_range) {
	  chunk_buffer_E_r(cur_ind) = 0.0;
	  cur_E_r_abs = 0.0;
	  num_truncations++;
	}
	if(cur_E_z_abs < max_field / chunkloop_data -> dynamic_range) {
	  chunk_buffer_E_z(cur_ind) = 0.0;
	  cur_E_z_abs = 0.0;
	}

	if((cur_E_r_abs > 0.0) && (cur_E_r_abs < min_abs_E_r)) {
	  min_abs_E_r = cur_E_r_abs;
	}
	if((cur_E_z_abs > 0.0) && (cur_E_z_abs < min_abs_E_z)) {
	  min_abs_E_z = cur_E_z_abs;
	}	
      }
      else {
	// The field has never gone above the threshold
	chunk_buffer_E_r(cur_ind) = 0.0;
	chunk_buffer_E_z(cur_ind) = 0.0;
	num_truncations++;
      }
    }

    std::cout << "truncated " << num_truncations << "/" << pts_r * pts_z << std::endl;
	
    chunkloop_data -> fstor.RegisterChunk(chunk_buffer_E_r, chunk_buffer_E_z, chunk_start_inds);        
  } // end chunkloop
  
} // end namespace meep
  
void CylindricalWeightingFieldCalculator::Calculate(std::filesystem::path outdir, std::filesystem::path mergedir, std::filesystem::path tmpdir) {

  // Prepare merge directory
  if(mergedir.empty()) {
    mergedir = outdir;
  }
  
  // TODO: this will get a lot easier once we can have all three components together in the same array
  std::filesystem::path mergedir_E_r = mergedir / "E_r";
  std::filesystem::path mergedir_E_z = mergedir / "E_z";
  
  if(meep::am_master()) {
    
    if(!std::filesystem::exists(mergedir)) {
      std::filesystem::create_directory(mergedir);
    }

    if(!std::filesystem::exists(mergedir_E_r)) {
      std::filesystem::create_directory(mergedir_E_r);
    }

    if(!std::filesystem::exists(mergedir_E_z)) {
      std::filesystem::create_directory(mergedir_E_z);
    }
  }
  
  // Make sure to have a valid temporary directory to store the results as we go along
  if(tmpdir.empty()) {
    char tmpdir_template[] = "/tmp/eisvogel.XXXXXX";
    tmpdir = std::filesystem::path(mkdtemp(tmpdir_template));
  }
  
  std::cout << "Using tmpdir = " << tmpdir << std::endl;
  std::cout << "Using mergedir = " << mergedir << std::endl;

  meep::all_wait();
  
  // This is only for cross-checking the geometry for now
  // f -> output_hdf5(meep::Dielectric, gv -> surroundings());

  IndexVector requested_chunk_size({400, 400, 400});
  std::shared_ptr<RZFieldStorage> fstor = std::make_shared<RZFieldStorage>(tmpdir, 10);

  int fstats_subsampling = 2;
  std::shared_ptr<FieldChunkStatisticsTracker2D> fstats = std::make_shared<FieldChunkStatisticsTracker2D>(fstats_subsampling);
  scalar_t dynamic_range = 50;
  scalar_t abs_min_field = 1e-20;
  ChunkloopData cld(0, *fstor, *fstats, dynamic_range, abs_min_field);

  // Setup simulation run
  m_f -> loop_in_chunks(meep::eisvogel_setup_chunkloop, static_cast<void*>(&cld), m_f -> total_volume());
  
  // Main simulation loop runs here  
  std::size_t stepcnt = 0;
  for(double cur_t = 0.0; cur_t <= m_t_end; cur_t += 0.1) {
  //for(double cur_t = 150.0; cur_t <= 152.0; cur_t += 0.1) {
  //for(double cur_t = 100.0; cur_t <= 102.0; cur_t += 0.13) {
  //for(double cur_t = 0.0; cur_t <= 1.0; cur_t += 0.13) {

    // Time-step the fields
    while (m_f -> time() < cur_t) {
      m_f -> step();
    }
    
    if(meep::am_master()) {
      std::cout << "Simulation time: " << m_f -> time() << std::endl;
    }
    
    cld.ind_t = stepcnt++;
    m_f -> loop_in_chunks(meep::eisvogel_saving_chunkloop, static_cast<void*>(&cld), m_f -> total_volume());

    std::cout << "LLL now on stepcnt = " << std::endl;
    std::cout << stepcnt << std::endl;
    
    if((stepcnt % 400) == 0) {
      std::cout << "BBB now merging chunks" << std::endl;
      fstor -> MergeChunks(0, 400);
    }
  }

  fstor -> MergeChunks(0, 400);
  
  // TODO: again, will get better once the three separate arrays are gone
  // TODO: for large weighting fields, will have to move chunks to the permanent location continuously throughout the calculation so as not to fill up local storage
  //       move them so that only the complete chunks (surviving after defragmentation) are moved that won't need to be accessed anymore
  std::cout << "moving output into merging location ...";
  std::filesystem::copy(tmpdir / "E_r", mergedir_E_r, std::filesystem::copy_options::recursive);
  std::filesystem::copy(tmpdir / "E_z", mergedir_E_z, std::filesystem::copy_options::recursive);
  std::cout << " done!" << std::endl;

  // Wait until everybody has finished copying
  meep::all_wait();

  std::cout << "==============================================" << std::endl;
  std::cout << "==============================================" << std::endl;
  std::cout << " All parallel things finished, now single-threaded defragmentation " << std::endl;
  std::cout << "==============================================" << std::endl;
  std::cout << "==============================================" << std::endl;
  
  if(meep::am_master()) {
    std::shared_ptr<CylindricalWeightingField> cwf = std::make_shared<CylindricalWeightingField>(mergedir, *m_start_coords, *m_end_coords);
    cwf -> MakeMetadataPersistent();    
    cwf -> RebuildChunks(requested_chunk_size);

    std::filesystem::copy(mergedir, outdir, std::filesystem::copy_options::recursive);
  }
}
