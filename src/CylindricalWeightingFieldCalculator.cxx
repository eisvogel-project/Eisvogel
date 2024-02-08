#include <iostream>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>
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
    m_subsampling(subsampling), m_largest_max(0.0), m_smallest_max(std::numeric_limits<scalar_t>::max()) { }

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

      if(cur_val > region_max_data(cur_subsampled_ind)) {
	
	// have a new local field maximum	
	region_max_data(cur_subsampled_ind) = cur_val;
	
	// also have a new global field maximum
	if(cur_val > m_largest_max) {
	  m_largest_max = cur_val;
	}
      }

      // do we have a new smallest non-zero maximum?
      if((cur_val > std::numeric_limits<scalar_t>::min()) && (cur_val < m_smallest_max)) {
	m_smallest_max = cur_val;
      }
    }    
  }

  scalar_t GetMaxLocal(KeyT& region_key, IndexVector& region_inds) {
    auto region_el = m_max_data.find(region_key);
    if(region_el == m_max_data.end()) {
      throw std::runtime_error("Error: requested region not found!");
    }    
    return region_el -> second(ToSubsampledInd(region_inds));
  }  

  scalar_t GetLargestMaxGlobal() {
    return m_largest_max;
  }

  scalar_t GetSmallestMaxGlobal() {
    return m_smallest_max;
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
  scalar_t m_largest_max;
  scalar_t m_smallest_max;
  
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

  ChunkloopData(std::size_t ind_t, RZFieldStorage& fstor, FieldChunkStatisticsTracker2D& fstats) :
    ind_t(ind_t), fstor(fstor), fstats(fstats) { }
  
  std::size_t ind_t;
  RZFieldStorage& fstor;
  FieldChunkStatisticsTracker2D& fstats;
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

    std::cout << "current global largest max: " << chunkloop_data -> fstats.GetLargestMaxGlobal() << std::endl;
    std::cout << "current global smallest max: " << chunkloop_data -> fstats.GetSmallestMaxGlobal() << std::endl;
    
    // auto to_keep = [](scalar_t value) -> bool {
    //   return std::fabs(value) > 1e-6;
    // };
    
    // SparseScalarField3D<scalar_t> sparse_chunk_buffer_E_r = SparseScalarField3D<scalar_t>::From(chunk_buffer_E_r, to_keep, 0.0);
    // SparseScalarField3D<scalar_t> sparse_chunk_buffer_E_z = SparseScalarField3D<scalar_t>::From(chunk_buffer_E_z, to_keep, 0.0);    
    // chunkloop_data -> fstor.RegisterChunk(sparse_chunk_buffer_E_r, sparse_chunk_buffer_E_z, chunk_start_inds);

    // chunkloop_data -> fstor.RegisterChunk(chunk_buffer_E_r, chunk_buffer_E_z, chunk_start_inds);
    
  } // end chunkloop
  
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

  int fstats_subsampling = 2;
  std::shared_ptr<FieldChunkStatisticsTracker2D> fstats = std::make_shared<FieldChunkStatisticsTracker2D>(fstats_subsampling);
  ChunkloopData cld(0, *fstor, *fstats);

  // Setup simulation run
  m_f -> loop_in_chunks(meep::eisvogel_setup_chunkloop, static_cast<void*>(&cld), m_f -> total_volume());
  
  // Main simulation loop runs here  
  std::size_t stepcnt = 0;
  for(double cur_t = 0.0; cur_t <= m_t_end; cur_t += 1) {

    // Time-step the fields
    while (m_f -> time() < cur_t) {
      m_f -> step();
    }
    
    if(meep::am_master()) {
      std::cout << "Simulation time: " << m_f -> time() << std::endl;
    }
    
    cld.ind_t = stepcnt++;
    std::cout << "entering saving chunkloop" << std::endl;
    m_f -> loop_in_chunks(meep::eisvogel_saving_chunkloop, static_cast<void*>(&cld), m_f -> total_volume());
    std::cout << "exit saving chunkloop" << std::endl;
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
