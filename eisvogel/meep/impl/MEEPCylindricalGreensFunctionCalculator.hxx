#include "Vector.hh"
#include "NDVecArray.hh"
#include "NDVecArrayOperations.hh"
#include "Symmetry.hh"
#include "GreensFunction.hh"
#include "ChunkUtils.hh"
#include <unordered_map>

template <typename RegionKeyT, std::size_t dims>
class FieldStatisticsTracker {

public:

  using region_t = NDVecArray<scalar_t, dims, 1>;
  using shape_t = typename region_t::shape_t;
  using ind_t = typename region_t::ind_t;
  
public:

  FieldStatisticsTracker(std::size_t downsampling) : m_downsampling(downsampling) {  }

  void AddRegion(const RegionKeyT& region_key, const shape_t& region_shape, scalar_t init_value = 0.0) {
    shape_t downsampled_shape = to_downsampled_shape(region_shape);
    m_max_data.emplace(std::make_pair(region_key, region_t(downsampled_shape, init_value)));

    std::cout << "job " << meep::my_rank() <<": for chunk with key = " << region_key << " now has " << m_max_data.size() << " chunks registered" << std::endl;
    // std::cout << "constructed max tracker with shape = " << subsampled_shape << std::endl;
    // std::cout << "for original chunk with shape = " << region_shape << std::endl;
  }

  void UpdateStatisticsForRegion(const RegionKeyT& region_key, const region_t& region_data) {

    // Update running field maximum
    assert(m_max_data.contains(region_key));
    region_t& region_max_data = m_max_data.at(region_key);

    auto max_updater = [&](const ind_t& ind){
      ind_t subsampled_ind = to_downsampled_ind(ind);
      scalar_t cur_val = region_data[ind][0];

      // have a new local field maximum; update
      if(cur_val > region_max_data[subsampled_ind][0]) {       
	region_max_data[subsampled_ind][0] = cur_val;
      }
    };    
    IteratorUtils::index_loop_over_array_elements(region_data, max_updater);
  }

  scalar_t GetMaxLocal(const RegionKeyT& region_key, const ind_t& ind) const {
    assert(m_max_data.contains(region_key));
    ind_t subsampled_ind = to_downsampled_ind(ind);
    return m_max_data.at(region_key)[subsampled_ind][0];
  }

private:
  
  shape_t to_downsampled_shape(const shape_t& shape) const {
    return (shape + m_downsampling - 1) / m_downsampling;
  }

  ind_t to_downsampled_ind(const ind_t& ind) const {
    return ind / m_downsampling;
  }
  
private:

  std::size_t m_downsampling;
  std::unordered_map<RegionKeyT, region_t> m_max_data;

};

// Data container to keep track of metadata information pertaining to a single spatial simulation chunk as used by MEEP.
template <typename SpatialSymmetryT>
struct SimulationChunkMetadata {

  using shape_t = Vector<std::size_t, SpatialSymmetryT::dims - 1>;
  using ind_t = Vector<std::size_t, SpatialSymmetryT::dims - 1>;

  SimulationChunkMetadata(const shape_t& chunk_shape, const ind_t& chunk_start_ind, const shape_t& storage_chunk_shape, const shape_t& storage_chunk_start_ind) :
    chunk_shape(chunk_shape), chunk_start_ind(chunk_start_ind), storage_chunk_shape(storage_chunk_shape), storage_chunk_start_ind(storage_chunk_start_ind) { }
  
  shape_t chunk_shape;  // Shape of this simulation chunk as it is used by MEEP
  ind_t chunk_start_ind;  // Start index of this simulation chunk

  shape_t storage_chunk_shape;  // Shape of the portion of this chunk that is to be saved
  ind_t storage_chunk_start_ind;  // Start index of the portion of the chunk that is to be saved
};

// Data container to pass into the MEEP callbacks defined below. Contains all the relevant information that needs to be passed
// back and forth between MEEP and Eisvogel.
template <typename SpatialSymmetryT>
struct ChunkloopData {

  using meep_chunk_ind_t = int;
  using sim_chunk_meta_t = SimulationChunkMetadata<SpatialSymmetryT>;
  using darr_t = typename SpatialSymmetryT::darr_t;
  using chunk_t = typename SpatialSymmetryT::darr_t::chunk_t;
  using spatial_chunk_t = typename SpatialSymmetryT::spatial_chunk_t;
  using view_t = typename chunk_t::view_t;  
  using fstats_t = FieldStatisticsTracker<meep_chunk_ind_t, SpatialSymmetryT::dims - 1>;

  ChunkloopData(std::size_t ind_time,
		spatial_chunk_t::shape_t calc_domain_shape, spatial_chunk_t::ind_t storage_domain_start_ind,
		spatial_chunk_t::shape_t storage_domain_shape,
		darr_t& fstor, fstats_t& fstats, scalar_t dynamic_range, std::size_t calc_chunk_size_t, std::size_t calc_chunk_size_linear,
		scalar_t abs_min_field, const chunk_t::ind_t& downsampling_factor,
		const chunk_t::shape_t& init_field_buffer_shape, const chunk_t::shape_t& init_field_buffer_processed_shape,
		const fstats_t::shape_t& init_stat_buffer_shape, const chunk_t::shape_t& init_field_chunk_buffer_shape) :
    ind_time(ind_time),
    calc_domain_shape(calc_domain_shape), storage_domain_start_ind(storage_domain_start_ind), storage_domain_shape(storage_domain_shape),
    fstor(fstor), fstats(fstats), dynamic_range(dynamic_range), calc_chunk_size_t(calc_chunk_size_t), calc_chunk_size_linear(calc_chunk_size_linear),
    abs_min_field(abs_min_field), downsampling_factor(downsampling_factor), field_buffer(init_field_buffer_shape),
    field_buffer_processed(init_field_buffer_processed_shape), field_absval_buffer(init_stat_buffer_shape),
    field_chunk_buffer(init_field_chunk_buffer_shape) { }
  
  std::size_t ind_time;  // Time index

  spatial_chunk_t::shape_t calc_domain_shape; // Shape of the region in which MEEP calculates the fields (always starts at zero)
  
  spatial_chunk_t::ind_t storage_domain_start_ind; // Start index of the domain in which the fields should be stored
  spatial_chunk_t::shape_t storage_domain_shape; // Shape of the domain in which the fields should be stored
  
  darr_t& fstor;  // Reference to field storage
  fstats_t& fstats;  // Reference to field statistics tracker
  
  scalar_t dynamic_range;  // Dynamic range of field to keep
  std::size_t calc_chunk_size_t;  // Chunk size along the time direction that is used during Green's function calculation
  std::size_t calc_chunk_size_linear;  // Chunk size along the spatial direction(s) used during Green's function calculation
  scalar_t abs_min_field;  // Abs. minimum field value to keep
  chunk_t::ind_t downsampling_factor;

  chunk_t field_buffer;  // Buffer to assemble field values from a single simulation chunk
  chunk_t field_buffer_processed;  // Buffer for field values with spatial downsampling applied
  fstats_t::region_t field_absval_buffer;  // Buffer for field magnitude (vector norm)
  chunk_t field_chunk_buffer;  // Buffer for storage chunks
  
  std::unordered_map<meep_chunk_ind_t, sim_chunk_meta_t> sim_chunk_meta;  // Metadata for all simulation chunks
};

// Local utilities
namespace {

  std::filesystem::path create_tmp_dir_in(std::filesystem::path parent_dir) {
    if(!std::filesystem::exists(parent_dir)) {
      std::filesystem::create_directory(parent_dir);
    }
    std::string tmpdir_template = parent_dir / "eisvogel.XXXXXX";    
    char* ret = mkdtemp(tmpdir_template.data());  (void)ret;
    return tmpdir_template;
  }

  // Ensures that `dir` is an empty directory, creating it if it does not yet exist,
  // or removing its contents if it does already exist
  void ensure_empty_directory(std::filesystem::path dir) {
    if(!std::filesystem::exists(dir)) {
      std::filesystem::create_directory(dir);
    }
    else {
      for(const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::filesystem::remove_all(entry.path());
      }
    }
  }

  template <class darr_t, class shape_t>
  bool verify_shape(std::filesystem::path dir, const shape_t& expected_shape) {
    darr_t darr(dir);
    std::cout << "GGGG found darr with shape = " << darr.GetShape() << ", expected " << expected_shape << std::endl;
    return(darr.GetShape() == expected_shape);
  }
}

using CylindricalChunkloopData = ChunkloopData<SpatialSymmetry::Cylindrical<scalar_t>>;

// Limit dynamic range
void limit_field_dynamic_range(int i_chunk, CylindricalChunkloopData::chunk_t& field_buffer, const CylindricalChunkloopData::fstats_t::region_t& field_absval_buffer,
			       const CylindricalChunkloopData::fstats_t& fstats, scalar_t abs_min_field, scalar_t dynamic_range) {

  using view_t = CylindricalChunkloopData::chunk_t::view_t; 

  std::size_t num_truncations = 0;

  auto truncator = [&](const TZRIndexVector& cur_ind, view_t cur_elem) {

    ZRIndexVector cur_slice_ind = cur_ind.zr();
    scalar_t max_abs_field = fstats.GetMaxLocal(i_chunk, cur_slice_ind);  // maximum field norm encountered so far at this locaiton

    if(max_abs_field > abs_min_field) {
      scalar_t cur_abs_field = field_absval_buffer[cur_slice_ind][0];  // current field norm at this location

      // truncate, we are above the required dynamic range
      if(cur_abs_field < max_abs_field / dynamic_range) {       
	cur_elem = 0;
	num_truncations++;
      }      
    }
    else {

      // The field has not yet gone above the `abs_min_field` threshold; truncate
      cur_elem = 0;
      num_truncations++;      
    }
  };  
  field_buffer.index_loop_over_elements(truncator);
}

// Small meep utility functions
namespace meep {

  ZRIndexVector get_cylindrical_chunk_start_ind(const ivec& is) {

    assert(is.z() >= 0);    
    assert(is.r() >= 0);

    return {
      (std::size_t)((is.z() - 1) / 2),
      (std::size_t)((is.r() - 1) / 2)
    };
  }

  ZRVector<std::size_t> get_cylindrical_chunk_shape(fields_chunk* fc, const ivec& is, const ivec& ie) {

    // determine rank (= number of spatial dimensions) of this simulation chunk
    std::size_t rank = number_of_directions(fc -> gv.dim);
    constexpr std::size_t required_rank = 2;
    assert(rank == required_rank);  // For a cylindrical geometry, we expect a 2d spatial simulation volume

    // determine the shape of this simulation chunk
    // shape = {shape_z, shape_r}
    std::size_t shape[required_rank] = {0};
    std::size_t index = 0;
    LOOP_OVER_DIRECTIONS(fc -> gv.dim, d) {
      std::size_t cur_len = std::max(0, (ie.in_direction(d) - is.in_direction(d)) / 2 + 1);
      shape[index++] = cur_len;
    }
    
    return {shape[0], shape[1]};    
  }
  
} // end namespace meep

// Callbacks to interface with MEEP
namespace meep {
  
  // This is called once at the beginning of the simulation run. It probes the chunk arrangement generated by meep
  // and sets up various things for later.
  void eisvogel_setup_cylindrical_storage_chunkloop(fields_chunk* fc, int ichunk, component, ivec is, ivec ie,
						    vec, vec, vec, vec, double, double,
						    ivec, std::complex<double>,
						    const symmetry&, int, void* cld) {
    
    CylindricalChunkloopData* chunkloop_data = static_cast<CylindricalChunkloopData*>(cld);

    ZRIndexVector chunk_start_ind = get_cylindrical_chunk_start_ind(is);
    ZRVector<std::size_t> chunk_shape = get_cylindrical_chunk_shape(fc, is, ie);

    // Nothing to be done for simulation chunks that don't overlap the region that is to be saved
    if(!ChunkUtils::overlaps(chunk_start_ind, chunk_shape,
			     chunkloop_data -> storage_domain_start_ind, chunkloop_data -> storage_domain_shape)) {

      std::cout << "------------" << std::endl;
      std::cout << "chunk #: " << ichunk << " does not overlap" << std::endl;
      std::cout << "------------" << std::endl;
      
      return;
    }
    
    // Determine the portion of this chunk that is to be stored
    ZRVector<std::size_t> storage_chunk_shape;
    ZRIndexVector storage_chunk_start_ind;
    ChunkUtils::get_chunk_overlap(chunk_start_ind, chunk_shape,
				  (ZRIndexVector)(chunkloop_data -> storage_domain_start_ind),
				  (ZRVector<std::size_t>)(chunkloop_data -> storage_domain_shape),
				  storage_chunk_start_ind, storage_chunk_shape);
    
    std::cout << "------------" << std::endl;
    std::cout << "chunk #: " << ichunk << std::endl;
    std::cout << "storage_domain_start_ind = " << chunkloop_data -> storage_domain_start_ind << std::endl;
    std::cout << "storage_domain_shape = " << chunkloop_data -> storage_domain_shape << std::endl;
    std::cout << "chunk_start_ind = " << chunk_start_ind << std::endl;
    std::cout << "chunk_shape = " << chunk_shape << std::endl;
    std::cout << "storage_chunk_start_ind = " << storage_chunk_start_ind << std::endl;
    std::cout << "storage_chunk_shape = " << storage_chunk_shape << std::endl;
    std::cout << "------------" << std::endl;
    
    CylindricalChunkloopData::sim_chunk_meta_t cur_meta(chunk_shape, chunk_start_ind, storage_chunk_shape, storage_chunk_start_ind);
    chunkloop_data -> sim_chunk_meta.emplace(std::make_pair(ichunk, cur_meta));
    
    // Setup field statistics tracker for this chunk
    chunkloop_data -> fstats.AddRegion(ichunk, storage_chunk_shape);
  }

  // This is called for every simulation timestep and performs the saving of the Green's function.
  void eisvogel_cylindrical_saving_chunkloop(fields_chunk* fc, int ichunk, component cgrid, ivec is, ivec ie,
					     vec, vec, vec, vec, double, double,
					     ivec shift, std::complex<double> shift_phase,
					     const symmetry& S, int sn, void* cld) {
    
    CylindricalChunkloopData* chunkloop_data = static_cast<CylindricalChunkloopData*>(cld);    

    // Nothing to be done for chunks that we don't need
    if(!(chunkloop_data -> sim_chunk_meta.contains(ichunk))) {
      return;
    }

    // std::cout << "------------" << std::endl;
    // std::cout << "chunk #: " << ichunk << std::endl;   
    
    // Make sure the buffers are of the correct size for this simulation chunk    
    ZRVector<std::size_t> spatial_storage_chunk_shape(chunkloop_data -> sim_chunk_meta.at(ichunk).storage_chunk_shape);
    TZRVector<std::size_t> field_slice_storage_shape(1, spatial_storage_chunk_shape);

    chunkloop_data -> field_absval_buffer.resize(spatial_storage_chunk_shape);
    chunkloop_data -> field_buffer.resize(field_slice_storage_shape);

    // Make sure that start index of this simulation chunk hasn't changed since the geometry was probed at the very beginning
    ZRIndexVector spatial_chunk_start_ind = get_cylindrical_chunk_start_ind(is);
    assert(spatial_chunk_start_ind == chunkloop_data -> sim_chunk_meta.at(ichunk).chunk_start_ind);

    ZRIndexVector spatial_storage_chunk_start_ind(chunkloop_data -> sim_chunk_meta.at(ichunk).storage_chunk_start_ind);
    TZRIndexVector storage_chunk_start_ind(chunkloop_data -> ind_time, spatial_storage_chunk_start_ind);

    // std::cout << "storage_chunk_start_ind = " << storage_chunk_start_ind << std::endl;
    
    // some preliminary setup
    vec rshift(shift * (0.5*fc->gv.inva));  // shift into unit cell for PBC geometries
    
    // prepare the list of field components to fetch at each grid point
    component components[] = {Ez, Er};
    chunkloop_field_components data(fc, cgrid, shift_phase, S, sn, 2, components);

    // loop over all grid points in chunk and fill the field buffer
    LOOP_OVER_IVECS(fc->gv, is, ie, idx) {

      // get grid indices and coordinates of parent point
      IVEC_LOOP_ILOC(fc->gv, iparent);  // grid indices
      IVEC_LOOP_LOC(fc->gv, rparent);   // cartesian coordinates
      
      // apply symmetry transform to get grid indices and coordinates of child point
      ivec ichild = S.transform(iparent, sn) + shift;
      vec rchild = S.transform(rparent, sn) + rshift;	    

      // Index of current simulation point
      ZRIndexVector cur_spatial_ind{
	(std::size_t)((ichild.z() - 1) / 2),
	(std::size_t)((ichild.r() - 1) / 2)
      };
      TZRIndexVector cur_ind(chunkloop_data -> ind_time, cur_spatial_ind);
      
      if(ChunkUtils::contains(spatial_storage_chunk_start_ind, spatial_storage_chunk_shape, cur_spatial_ind)) {

	// fetch field components at child point ...
	data.update_values(idx);
	double E_z_val = data.values[0].real();
	double E_r_val = data.values[1].real();      
	double E_abs_val = std::sqrt(E_z_val * E_z_val + E_r_val * E_r_val);
	
	// ... and store them
	CylindricalChunkloopData::view_t field_elem = chunkloop_data -> field_buffer[cur_ind - storage_chunk_start_ind];
	field_elem[0] = (scalar_t)E_r_val;
	field_elem[1] = (scalar_t)E_z_val;
	chunkloop_data -> field_absval_buffer[cur_spatial_ind - spatial_storage_chunk_start_ind] = (scalar_t)E_abs_val;	
      }      
    }

    // Update statistics tracker
    chunkloop_data -> fstats.UpdateStatisticsForRegion(ichunk, chunkloop_data -> field_absval_buffer);

    // Truncate small field values to keep a certain maximum dynamic range
    limit_field_dynamic_range(ichunk, chunkloop_data -> field_buffer, chunkloop_data -> field_absval_buffer, chunkloop_data -> fstats,
			      chunkloop_data -> abs_min_field, chunkloop_data -> dynamic_range);

    // Downsample the field buffer before it is registered in the output storage
    // Start index of this storage chunk relative to the full storage domain
    TZRIndexVector storage_chunk_rel_start_ind(chunkloop_data -> ind_time,
					       spatial_storage_chunk_start_ind - chunkloop_data -> storage_domain_start_ind);

    // std::cout << "storage_chunk_rel_start_ind = " << storage_chunk_rel_start_ind << std::endl;
    // std::cout << "downsampling_factor = " << chunkloop_data -> downsampling_factor << std::endl;
    
    TZRIndexVector processed_chunk_start_ind(0);
    Downsampling::downsample(chunkloop_data -> field_buffer, storage_chunk_rel_start_ind, chunkloop_data -> downsampling_factor,
			     chunkloop_data -> field_buffer_processed, processed_chunk_start_ind);
    assert((chunkloop_data -> field_buffer_processed).GetShape() == Downsampling::get_downsampled_shape(storage_chunk_rel_start_ind,
													field_slice_storage_shape,
													chunkloop_data -> downsampling_factor));

    // std::cout << "processed_chunk_start_ind = " << processed_chunk_start_ind << std::endl;
    
    // The downsampling has resulted in an empty chunk, nothing further to do
    if((chunkloop_data -> field_buffer_processed).GetNumberElements() == 0) {
      return;
    }

    // std::cout << "chunk shape before downsampling = " << chunkloop_data -> field_buffer.GetShape() << ", after downsampling = " <<
    //   chunkloop_data -> field_buffer_processed.GetShape() << std::endl;
    
    using shape_t = CylindricalChunkloopData::chunk_t::shape_t;
    shape_t requested_output_chunk_size(chunkloop_data -> calc_chunk_size_linear);

    // Register the processed buffer
    if(chunkloop_data -> ind_time % chunkloop_data -> calc_chunk_size_t == 0) {
      
      // Register this slice as the beginning of a new chunk ...
      TZRIndexVector start(0);
      auto simulation_chunk_storer = [&](const TZRIndexVector& output_chunk_start, const TZRIndexVector& output_chunk_end) {
	chunkloop_data -> field_chunk_buffer.resize(output_chunk_end - output_chunk_start);
	chunkloop_data -> field_chunk_buffer.fill_from(chunkloop_data -> field_buffer_processed,
						       output_chunk_start,
						       output_chunk_end,
						       start);
	chunkloop_data -> fstor.RegisterChunk(chunkloop_data -> field_chunk_buffer, processed_chunk_start_ind + output_chunk_start);
      };
      IteratorUtils::index_loop_over_chunks(start, chunkloop_data -> field_buffer_processed.GetShape(), requested_output_chunk_size,
					    simulation_chunk_storer);
      
      // std::cout << "registering new chunk at start_ind = " << chunk_start_ind << " with shape = " << chunkloop_data -> field_buffer.GetShape() << std::endl;      
    }
    else {

      // ... or append it to an already-existing chunk along the outermost (time) direction
      TZRIndexVector start(0);
      auto simulation_chunk_appender = [&](const TZRIndexVector& output_chunk_start, const TZRIndexVector& output_chunk_end) {
	chunkloop_data -> field_chunk_buffer.resize(output_chunk_end - output_chunk_start);
	chunkloop_data -> field_chunk_buffer.fill_from(chunkloop_data -> field_buffer_processed,
						       output_chunk_start,
						       output_chunk_end,
						       start);
	chunkloop_data -> fstor.AppendSlice<0>(processed_chunk_start_ind + output_chunk_start, chunkloop_data -> field_chunk_buffer);
      };
      IteratorUtils::index_loop_over_chunks(start, chunkloop_data -> field_buffer_processed.GetShape(), requested_output_chunk_size,
					    simulation_chunk_appender);
      
      // std::cout << "appending slice at start_ind = " << chunk_start_ind << " with shape = " << chunkloop_data -> field_buffer.GetShape() << std::endl;     
    }

    // std::cout << "------------" << std::endl;
  }  
} // end namespace meep

namespace GreensFunctionCalculator::MEEP {

  struct CylindricalGreensFunctionCalculator::MPIChunkCalculationResult {
    
    MPIChunkCalculationResult(CylinderRegion& region_stored, ZRVector<std::size_t>& padding_pre, ZRVector<std::size_t>& padding_post) :
      region_stored(region_stored), padding_pre(padding_pre), padding_post(padding_post) { }
    
    CylinderRegion region_stored;  // region of the cylindrical simulation volume that has been stored to disk
    ZRVector<std::size_t> padding_pre;  // spatial padding on the `bottom-left' of the calculated domain (number of elements)
    ZRVector<std::size_t> padding_post;  // spatial padding on the `top-right' of the calculated domain (number of elements)
  };
  
  CylindricalGreensFunctionCalculator::MPIChunkCalculationResult
  CylindricalGreensFunctionCalculator::calculate_mpi_chunk(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::size_t padding_requested,
							   double courant_factor, double resolution, double timestep, std::size_t calc_chunk_size_t,
							   std::size_t calc_chunk_size_linear, double pml_width, std::size_t downsampling_on_disk,
							   scalar_t dynamic_range, scalar_t abs_min_field) {
    
    // Make sure that MEEP will not silently round the domain size
    if(!MathUtils::is_integer(m_geom.GetRSize() * resolution)) {
      throw std::runtime_error("Error: chosen radial size of geometry would introduce rounding.");
    }

    if(!MathUtils::is_integer(m_geom.GetZSize() * resolution)) {
      throw std::runtime_error("Error: chosen axial size of geometry would introduce rounding.");
    }
    
    std::shared_ptr<meep::grid_volume> meep_gv = std::make_shared<meep::grid_volume>(meep::volcyl(m_geom.GetRMax(), m_geom.GetZMax() - m_geom.GetZMin(), resolution));
    std::shared_ptr<meep::structure> meep_s = std::make_shared<meep::structure>(*meep_gv, m_geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
    std::shared_ptr<meep::fields> meep_f = std::make_shared<meep::fields>(meep_s.get());

    // Get the absolute value of the Green's function correct
    scalar_t magic_number = 0.11f;

    // A "point source" in MEEP has a size of 1 voxel, i.e. 1 / resolution in MEEP units,
    // but later in the integration code we assume ds = 1 as the size of the infinitesimal dipole
    scalar_t feed_current_scale_factor = resolution * magic_number;
    m_antenna.SetScaleFactor(feed_current_scale_factor);
    
    m_antenna.AddToGeometry(*meep_f, m_geom);
    
    // Some typedefs
    constexpr std::size_t dims = SpatialSymmetry::Cylindrical<scalar_t>::dims;
    // using chunk_t = typename SpatialSymmetry::Cylindrical<scalar_t>::chunk_t;
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
    
    // Working directory in process-local scratch area
    std::filesystem::path local_workdir = create_tmp_dir_in(local_scratchdir);  
    
    std::size_t cache_depth = 0;   // all chunks are created one time-slice at a time, direcly request them to be streamed to disk
    TZRVector<std::size_t> init_cache_el_shape(1);
    TZRVector<std::size_t> streamer_chunk_shape(stor::INFTY); streamer_chunk_shape[0] = 1; // serialize one time slice at a time
    darr_t darr(local_workdir, cache_depth, init_cache_el_shape, streamer_chunk_shape);
    
    // Set up tracker to follow field statistics such as the maximum field strength
    // This corresponds to a constant-time slice, which has one fewer dimensions
    // (`fstats_downsampling` is unrelated to `downsampling_on_disk`; the former
    // affects only the field statistics tracker, the latter affects only the chunk
    // that is stored)
    std::size_t fstats_downsampling = 2;  
    FieldStatisticsTracker<meep_chunk_ind_t, dims - 1> fstats(fstats_downsampling);
        
    // Figure out which index range gives us the region that we are asked to store
    ZRVector<std::size_t> calc_domain_shape{(std::size_t)(meep_gv -> nz()), (std::size_t)(meep_gv -> nr())};
    ZRIndexVector storage_domain_start_ind{
      std::size_t((m_request_to_store.GetZMin() - m_geom.GetZMin()) / m_geom.GetZSize() * calc_domain_shape.z()),
      std::size_t((m_request_to_store.GetRMin() - m_geom.GetRMin()) / m_geom.GetRSize() * calc_domain_shape.r())
    };
    ZRVector<std::size_t> storage_domain_shape{
      std::size_t(m_request_to_store.GetZSize() / m_geom.GetZSize() * calc_domain_shape.z()),
      std::size_t(m_request_to_store.GetRSize() / m_geom.GetRSize() * calc_domain_shape.r()),
    };

    // Communicate back which region we are actually storing now
    assert(calc_domain_shape.z() == m_geom.GetZSize() * resolution);
    assert(calc_domain_shape.r() == m_geom.GetRSize() * resolution);
    
    ZRIndexVector storage_domain_end_ind = storage_domain_start_ind + storage_domain_shape;

    // std::cout << "HHHH region_request: r_min = " << m_request_to_store.GetRMin() << std::endl;
    // std::cout << "HHHH region_request: r_max = " << m_request_to_store.GetRMax() << std::endl;
    // std::cout << "HHHH region_request: z_min = " << m_request_to_store.GetZMin() << std::endl;
    // std::cout << "HHHH region_request: z_max = " << m_request_to_store.GetZMax() << std::endl;
    
    // std::cout << "HHHH calc_domain_shape = " << calc_domain_shape << std::endl;
    // std::cout << "HHHH storage_domain_shape = " << storage_domain_shape << std::endl;
    // std::cout << "HHHH storage_domain_start_ind = " << storage_domain_start_ind << std::endl;
    // std::cout << "HHHH storage_domain_end_ind = " << storage_domain_end_ind << std::endl;    
    // std::cout << "setting region" << std::endl;

    CylinderRegion region_stored((scalar_t)(storage_domain_start_ind.r()) / resolution + m_geom.GetRMin(),
				 (scalar_t)(storage_domain_end_ind.r()) / resolution + m_geom.GetRMin(),
				 (scalar_t)(storage_domain_start_ind.z()) / resolution + m_geom.GetZMin(),
				 (scalar_t)(storage_domain_end_ind.z()) / resolution + m_geom.GetZMin());

    // Attempt to extend the storage domain by adding the requested padding around it
    // Note: `padding_pre` and `padding_post` here correspond to the padding that will exist *after* the downsampling has been applied
    ZRVector<std::size_t> padding_pre = VectorUtils::min(storage_domain_start_ind / downsampling_on_disk, padding_requested);
    ZRVector<std::size_t> padding_post = VectorUtils::min((calc_domain_shape - storage_domain_end_ind) / downsampling_on_disk, padding_requested);
    assert(storage_domain_start_ind >= padding_pre * downsampling_on_disk);
    assert(calc_domain_shape >= storage_domain_end_ind + padding_post * downsampling_on_disk);

    // std::cout << "HHHH padding_pre = " << padding_pre << std::endl;
    // std::cout << "HHHH padding_post = " << padding_post << std::endl;
    
    // Prepare data container to pass to all MEEP callbacks
    TZRVector<std::size_t> downsampling_factor(downsampling_on_disk); downsampling_factor.t() = 1;  // downsample only along the spatial directions
    TZRVector<std::size_t> init_field_buffer_shape(1);
    TZRVector<std::size_t> init_field_chunk_buffer_shape(1);
    ZRVector<std::size_t> init_field_absval_buffer_shape(1);

    // Now extend the storage domain such that, after downsampling, the correct padding is achieved
    ZRIndexVector storage_domain_padded_start_ind = storage_domain_start_ind - padding_pre * downsampling_on_disk;
    ZRVector<std::size_t> storage_domain_padded_shape = storage_domain_shape + (padding_pre + padding_post) * downsampling_on_disk;

    // std::cout << "HHH expected shape with padding after downsampling: " << Downsampling::get_downsampled_shape(ZRVector<std::size_t>(0), storage_domain_padded_shape, ZRVector<std::size_t>(downsampling_on_disk)) << std::endl;
    // std::cout << "HHH expected shape without padding after downsampling: " << Downsampling::get_downsampled_shape(ZRVector<std::size_t>(0), storage_domain_shape, ZRVector<std::size_t>(downsampling_on_disk)) << std::endl;
    
    assert(Downsampling::get_downsampled_shape(ZRVector<std::size_t>(0), storage_domain_padded_shape, ZRVector<std::size_t>(downsampling_on_disk)) ==
	   Downsampling::get_downsampled_shape(ZRVector<std::size_t>(0), storage_domain_shape, ZRVector<std::size_t>(downsampling_on_disk))
	   + padding_pre + padding_post);
    
    // Package up all the information
    CylindricalChunkloopData cld(0, calc_domain_shape, storage_domain_padded_start_ind, storage_domain_padded_shape,
				 darr, fstats, dynamic_range, calc_chunk_size_t, calc_chunk_size_linear, abs_min_field, downsampling_factor,
				 init_field_buffer_shape, init_field_buffer_shape,
				 init_field_absval_buffer_shape, init_field_chunk_buffer_shape);

    // std::cout << "HHHH storage_domain_padded_start_ind = " << storage_domain_padded_start_ind << std::endl;
    // std::cout << "HHHH storage_domain_padded_shape = " << storage_domain_padded_shape << std::endl;
    
    // std::cout << "HHHH region_stored: r_min = " << region_stored.GetRMin() << std::endl;
    // std::cout << "HHHH region_stored: r_max = " << region_stored.GetRMax() << std::endl;
    // std::cout << "HHHH region_stored: z_min = " << region_stored.GetZMin() << std::endl;
    // std::cout << "HHHH region_stored: z_max = " << region_stored.GetZMax() << std::endl;
    
    // Setup simulation run
    meep_f -> loop_in_chunks(meep::eisvogel_setup_cylindrical_storage_chunkloop, static_cast<void*>(&cld), meep_f -> total_volume());
    
    // Main simulation loop
    std::size_t stepcnt = 0;
    for(double cur_t = 0.0; cur_t <= m_t_end; cur_t += timestep) {
      
      // Time-step the fields
      while (meep_f -> time() < cur_t) {
	meep_f -> step();
      }
      
      if(meep::am_master()) {
	std::cout << "Simulation time: " << meep_f -> time() << std::endl;
      }
      
      // Store the Green's function at the present timestep
      cld.ind_time = stepcnt++;
      meep_f -> loop_in_chunks(meep::eisvogel_cylindrical_saving_chunkloop, static_cast<void*>(&cld), meep_f -> total_volume());
    }
    
    // Swap axes 0 and 2 to go from TZR to RZT
    darr.SwapAxes<0, 2>();
    
    // Move Green's function to the requested output location
    darr.Move(outdir);

    // Return information describing the result of this calculation
    return MPIChunkCalculationResult(region_stored, padding_pre, padding_post);
  }

  // No dedicated storage region is passed, store everything
  CylindricalGreensFunctionCalculator::CylindricalGreensFunctionCalculator(CylinderGeometry& geom, Antenna& antenna, scalar_t t_end) :
    CylindricalGreensFunctionCalculator(geom, CylinderRegion((scalar_t)(0.0), geom.GetRMax(), geom.GetZMin(), geom.GetZMax()), antenna, t_end) { }

  CylindricalGreensFunctionCalculator::CylindricalGreensFunctionCalculator(CylinderGeometry& geom, const CylinderRegion& request_to_store,
									   Antenna& antenna, scalar_t t_end) :
    m_t_end(t_end), m_geom(geom), m_request_to_store(request_to_store), m_antenna(antenna) {

    // Sanity checks
    if(!m_geom.contains(m_request_to_store)) {
      throw std::runtime_error("Error: requested storage region contains undefined domain!");
    }

    if(m_request_to_store.IsEmpty()) {
      throw std::runtime_error("Error: requested storage region is empty!");
    }
  }
  
  void CylindricalGreensFunctionCalculator::merge_mpi_chunks(std::filesystem::path outdir, const std::vector<std::filesystem::path>& indirs) {
    
    // This must only be run by one process
    assert(meep::am_master());
    
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
    
    ensure_empty_directory(outdir);
    
    std::size_t cache_depth = 0;
    TZRVector<std::size_t> init_cache_el_shape(1);
    TZRVector<std::size_t> streamer_chunk_shape(stor::INFTY); streamer_chunk_shape[0] = 1; // serialize one time slice at a time
    darr_t darr(outdir, cache_depth, init_cache_el_shape, streamer_chunk_shape); // no big cache needed here, this is just for merging
    
    for(const std::filesystem::path& indir : indirs) {
      darr.Import(indir, true);
      std::filesystem::remove_all(indir);
    }
  }
  
  // Tries to come up with a good `partition_size` so that `partition_size` is an integer multiple of `chunk_size` and it takes
  // approximately `number_partitions_requested` nonoverlapping partitions to cover a region with `shape`.
  template <std::size_t dims>
  Vector<std::size_t, dims> get_partition_size(const Vector<std::size_t, dims>& shape, const Vector<std::size_t, dims>& chunk_size,
					       std::size_t number_partitions_requested) {
    
    std::size_t longest_axis = std::distance(shape.begin(), std::max_element(std::execution::unseq, shape.begin(), shape.end()));
    
    auto number_partitions_actual = [&](const Vector<std::size_t, dims>& pars_per_dim) -> std::size_t {
      return std::reduce(std::execution::unseq, pars_per_dim.begin(), pars_per_dim.end(), 1, std::multiplies<int>());
    };  
    
    // Number of (partial) chunks along each direction
    Vector<std::size_t, dims> chunks_per_dim = VectorUtils::ceil_nonneg(shape.template as_type<scalar_t>() / chunk_size.template as_type<scalar_t>());
    
    // First guess for the number of partitions per dimension
    Vector<std::size_t, dims> pars_per_dim((std::size_t)std::pow((scalar_t)number_partitions_requested, 1.0 / dims));
    
    while(number_partitions_actual(pars_per_dim) < number_partitions_requested) {
    pars_per_dim[longest_axis]++;
    }  
    
    Vector<std::size_t, dims> par_size_in_chunks = VectorUtils::ceil_nonneg(chunks_per_dim.template as_type<scalar_t>() / pars_per_dim.template as_type<scalar_t>());
    Vector<std::size_t, dims> partition_size = par_size_in_chunks * chunk_size;
    
    return partition_size;
  }
  
  void CylindricalGreensFunctionCalculator::rechunk_mpi(std::filesystem::path outdir, std::filesystem::path indir, std::filesystem::path global_scratch_dir,
							int cur_mpi_id, int number_mpi_jobs, const RZTVector<std::size_t>& requested_chunk_size,
							std::size_t overlap, const RZTVector<std::size_t>& padding_pre, const RZTVector<std::size_t>& padding_post,
							std::size_t cache_depth) {
    
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
    
    RZTVector<std::size_t> init_cache_el_shape(1);
    RZTVector<std::size_t> streamer_chunk_shape(stor::INFTY); streamer_chunk_shape[0] = 1; // serialize one radial slice at a time
    darr_t darr(indir, cache_depth, init_cache_el_shape, streamer_chunk_shape);
    
    std::cout << "loading array from " << indir << std::endl;
    std::cout << "now rechunking; have input array with shape " << darr.GetShape() << std::endl;
    std::cout << "this is rechunking job " << cur_mpi_id << " / " << number_mpi_jobs << " jobs" << std::endl;
    std::cout << "rechunking onto chunk_size = " << requested_chunk_size << " and overlap = " << overlap << std::endl;

    // Build the domain that will be rechunked, considering any padding that has been added before
    RZTIndexVector rechunk_start(padding_pre);
    RZTIndexVector rechunk_shape(darr.GetShape() - padding_pre - padding_post);
    RZTIndexVector rechunk_end = rechunk_start + rechunk_shape;

    std::cout << "after rechunking and accounting for padding will have array with shape " << rechunk_shape << ", starting at " << rechunk_start << std::endl;
    
    // Determine nonoverlapping rectangular partitions that are handled by each job; try to make them of similar size for optimal work sharing
    RZTVector<std::size_t> partition_size = get_partition_size(rechunk_shape, requested_chunk_size, number_mpi_jobs);
    
    auto rechunk_dir = [](int rechunk_num) -> std::filesystem::path {
      return std::format("rechunk_job_{}", rechunk_num);
    };
    
    int num_rechunking = 0;
    std::vector<std::filesystem::path> rechunking_dirs;
    auto rechunker = [&](const RZTIndexVector& partition_start, const RZTIndexVector& partition_end) {
      
      std::filesystem::path cur_rechunking_dir = global_scratch_dir / rechunk_dir(num_rechunking);
      rechunking_dirs.push_back(cur_rechunking_dir);
      
      // Each job only performs those rechunkings assigned to it
      if(num_rechunking % number_mpi_jobs == cur_mpi_id) {
	
	ensure_empty_directory(cur_rechunking_dir);
	
	std::cout << "job " << cur_mpi_id << " rechunking " << partition_start << " -> " << partition_end << " into " << cur_rechunking_dir << std::endl;
	darr.RebuildChunksPartial(partition_start, partition_end,
				  partition_start - rechunk_start, // The rechunked output array should be zero-indexed again
				  requested_chunk_size, cur_rechunking_dir, overlap, SpatialSymmetry::Cylindrical<scalar_t>::boundary_evaluator,
				  ChunkHints::ENABLE_OPT);
	std::cout << "job " << cur_mpi_id << " done" << std::endl;
      }
      
      num_rechunking++;
    };
    IteratorUtils::index_loop_over_chunks(rechunk_start, rechunk_end, partition_size, rechunker);
    
    // After all jobs are finished with rechunking their respective regions, merge all outputs into `outdir`
    meep::all_wait(); 
    
    if(meep::am_master()) {
      merge_mpi_chunks(outdir, rechunking_dirs);

      // Verify that the shape is what is expected
      assert(verify_shape<darr_t>(outdir, rechunk_shape));
    }
  }
  
  void CylindricalGreensFunctionCalculator::Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
						      double courant_factor, double resolution, double timestep, double pml_width, std::size_t downsampling_on_disk,
						      scalar_t dynamic_range, scalar_t abs_min_field, std::size_t chunk_overlap, std::size_t chunk_size_linear,
						      std::size_t rechunk_cache_depth, std::size_t calc_chunk_size_t, std::size_t calc_chunk_size_linear) {
    
    int number_mpi_jobs = meep::count_processors();
    int job_mpi_rank = meep::my_rank();
    
    // Prepare an empty working directory in global scratch area
    std::filesystem::path global_workdir = global_scratchdir / "eisvogel.global";
    if(meep::am_master()) {
      if(std::filesystem::exists(global_workdir)) {
	std::filesystem::remove_all(global_workdir);
      }
      std::filesystem::create_directory(global_workdir);
    }
    
    // Output directory
    if(meep::am_master()) {
      if(std::filesystem::exists(outdir) && !std::filesystem::is_empty(outdir)) {
	throw std::runtime_error("Error: cowardly refusing to overwrite an already-existing Green's function!");
      }
      std::filesystem::create_directory(outdir);
    }
    
    std::cout << "using global_workdir = " << global_workdir << std::endl;
    
    meep::all_wait();
    
    auto calc_job_dir = [](int job_mpi_rank) -> std::filesystem::path {
      return std::format("calc_job_{}", job_mpi_rank);
    };
    
    // First stage: each MPI process calculates and stores its part of the Green's function in `job_outdir`
    std::filesystem::path job_outdir = global_workdir / calc_job_dir(job_mpi_rank);
    std::size_t padding_requested = chunk_overlap; // number of additional elements that should be foreseen as padding around the calculated array
    MPIChunkCalculationResult mpi_calc_res = calculate_mpi_chunk(job_outdir, local_scratchdir, padding_requested,
								 courant_factor, resolution, timestep, calc_chunk_size_t, calc_chunk_size_linear,
								 pml_width, downsampling_on_disk, dynamic_range, abs_min_field);
    
    meep::all_wait();   // Wait for all processes to have finished providing their output
    
    // Second stage: create a new distributed array and import all the chunks. This is now a complete Green's function!
    std::filesystem::path mergedir = global_workdir / "merge";      
    if(meep::am_master()) {
      
      // Assemble output directories from all jobs
      std::vector<std::filesystem::path> to_merge;    
      for(int job_id = 0; job_id < number_mpi_jobs; job_id++) {
	std::filesystem::path cur_job_outdir = global_workdir / calc_job_dir(job_id);
	to_merge.push_back(cur_job_outdir);
      }
      
      // Call the merger
      merge_mpi_chunks(mergedir, to_merge);
    }
    
    meep::all_wait();
    
    // Third stage: rechunk the complete array and introduce overlap between neighbouring chunks, if requested
    RZTVector<std::size_t> requested_chunk_size(chunk_size_linear);
    RZTVector<std::size_t> padding_pre(mpi_calc_res.padding_pre.r(), mpi_calc_res.padding_pre.z(), 0u);
    RZTVector<std::size_t> padding_post(mpi_calc_res.padding_post.r(), mpi_calc_res.padding_post.z(), 0u);
    rechunk_mpi(outdir, mergedir, global_workdir, job_mpi_rank, number_mpi_jobs, requested_chunk_size,
		chunk_overlap, padding_pre, padding_post, rechunk_cache_depth);
    
    meep::all_wait();
    
    // Clean up our temporary files
    if(meep::am_master()) {
      std::filesystem::remove_all(global_workdir);
      
      // Fourth stage: create the actual Green's function object
      using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
      darr_t darr(outdir);
      RZTVector<std::size_t> darr_shape = darr.GetShape();

      RZTCoordVector storage_start_coords(mpi_calc_res.region_stored.GetRMin(), mpi_calc_res.region_stored.GetZMin(), (scalar_t)0.0);
      RZTCoordVector storage_end_coords(mpi_calc_res.region_stored.GetRMax(), mpi_calc_res.region_stored.GetZMax(), m_t_end);
            
      std::cout << "Green's function with storage_start = " << storage_start_coords << " and storage_end = " << storage_end_coords << std::endl;
      
      RZTCoordVector step_size = (storage_end_coords - storage_start_coords) / (darr_shape.template as_type<scalar_t>() - 1);
      CylindricalGreensFunction(storage_start_coords, storage_end_coords, step_size, std::move(darr));
    }
  }
}
