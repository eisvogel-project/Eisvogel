#include "Vector.hh"
#include "NDVecArray.hh"
#include "Symmetry.hh"
#include "GreensFunction.hh"
#include <unordered_map>

// Utilities for downsampling
namespace {

  template <typename ShapeT>
  ShapeT to_downsampled_shape(const ShapeT& shape, std::size_t downsampling) {
    return (shape + downsampling - 1) / downsampling;
  }

  template <typename IndT>
  IndT to_downsampled_ind(const IndT& ind, std::size_t downsampling) {
    return ind / downsampling;
  }  
}

template <typename RegionKeyT, std::size_t dims>
class FieldStatisticsTracker {

public:

  using region_t = NDVecArray<scalar_t, dims, 1>;
  using shape_t = typename region_t::shape_t;
  using ind_t = typename region_t::ind_t;
  
public:

  FieldStatisticsTracker(std::size_t downsampling) : m_downsampling(downsampling) {  }

  void AddRegion(const RegionKeyT& region_key, const shape_t& region_shape, scalar_t init_value = 0.0) {
    shape_t subsampled_shape = to_downsampled_shape(region_shape, m_downsampling);
    m_max_data.emplace(std::make_pair(region_key, region_t(subsampled_shape, init_value)));

    std::cout << "job " << meep::my_rank() <<": for chunk with key = " << region_key << " now has " << m_max_data.size() << " chunks registered" << std::endl;
    // std::cout << "constructed max tracker with shape = " << subsampled_shape << std::endl;
    // std::cout << "for original chunk with shape = " << region_shape << std::endl;
  }

  void UpdateStatisticsForRegion(const RegionKeyT& region_key, const region_t& region_data) {

    // Update running field maximum
    assert(m_max_data.contains(region_key));
    region_t& region_max_data = m_max_data.at(region_key);

    auto max_updater = [&](const ind_t& ind){
      ind_t subsampled_ind = to_downsampled_ind(ind, m_downsampling);
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
    ind_t subsampled_ind = to_downsampled_ind(ind, m_downsampling);
    return m_max_data.at(region_key)[subsampled_ind][0];
  }
  
private:

  std::size_t m_downsampling;
  std::unordered_map<RegionKeyT, region_t> m_max_data;

};

// Data container to keep track of metadata information pertaining to a single spatial simulation chunk.
template <typename SpatialSymmetryT>
struct SimulationChunkMetadata {

  using shape_t = Vector<std::size_t, SpatialSymmetryT::dims - 1>;

  SimulationChunkMetadata(const shape_t& chunk_shape, const shape_t& downsampled_chunk_shape) : chunk_shape(chunk_shape), downsampled_chunk_shape(downsampled_chunk_shape) { }
  
  shape_t chunk_shape;  // Shape of this simulation chunk as it is used by MEEP
  shape_t downsampled_chunk_shape;  // Shape of this simulation chunk after downsampling is applied
};

// Data container to pass into the MEEP callbacks defined below. Contains all the relevant information that needs to be passed
// back and forth between MEEP and Eisvogel.
template <typename SpatialSymmetryT>
struct ChunkloopData {

  using meep_chunk_ind_t = int;
  using sim_chunk_meta_t = SimulationChunkMetadata<SpatialSymmetryT>;
  using darr_t = typename SpatialSymmetryT::darr_t;
  using chunk_t = typename SpatialSymmetryT::darr_t::chunk_t;
  using view_t = typename chunk_t::view_t;  
  using fstats_t = FieldStatisticsTracker<meep_chunk_ind_t, SpatialSymmetryT::dims - 1>;

  ChunkloopData(std::size_t ind_time, darr_t& fstor, fstats_t& fstats, scalar_t dynamic_range, scalar_t abs_min_field, std::size_t downsampling_on_disk,
		const chunk_t::shape_t& init_field_buffer_shape, const chunk_t::shape_t& init_field_buffer_processed_shape,
		const fstats_t::shape_t& init_stat_buffer_shape, const chunk_t::shape_t& init_field_chunk_buffer_shape) :
    ind_time(ind_time), fstor(fstor), fstats(fstats), dynamic_range(dynamic_range), abs_min_field(abs_min_field), downsampling_on_disk(downsampling_on_disk),
    field_buffer(init_field_buffer_shape), field_buffer_processed(init_field_buffer_processed_shape),
    field_absval_buffer(init_stat_buffer_shape), field_chunk_buffer(init_field_chunk_buffer_shape) { }
  
  std::size_t ind_time;  // Time index
  darr_t& fstor;  // Reference to field storage
  fstats_t& fstats;  // Reference to field statistics tracker
  
  scalar_t dynamic_range;  // Dynamic range of field to keep
  scalar_t abs_min_field;  // Abs. minimum field value to keep
  std::size_t downsampling_on_disk;

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

// Applies downsampling to `field_buffer` along all axes
void downsample_field(const CylindricalChunkloopData::chunk_t& field_buffer, CylindricalChunkloopData::chunk_t& field_buffer_downsampled, const std::size_t downsampling) {
  field_buffer_downsampled.resize(to_downsampled_shape(field_buffer.GetShape(), downsampling));

  using ind_t = typename CylindricalChunkloopData::chunk_t::ind_t;

  auto downsampler = [&](const ind_t& ind_to_keep, const ind_t&) {
    ind_t downsampled_ind = to_downsampled_ind(ind_to_keep, downsampling);
    field_buffer_downsampled[downsampled_ind] = field_buffer[ind_to_keep];
  };  
  ind_t start(0);
  ind_t downsampling_chunk_shape(downsampling); // shape of region for which only one sample (the first one) needs to be kept
  IteratorUtils::index_loop_over_chunks(start, field_buffer.GetShape(), downsampling_chunk_shape, downsampler);
}

// Callbacks to interface with MEEP
namespace meep {

  // This is called once at the very beginning of the simulation run. It probes the chunk arrangement generated by meep
  // and sets up various things for later.
  void eisvogel_setup_chunkloop(fields_chunk* fc, int ichunk, component, ivec is, ivec ie,
				vec, vec, vec, vec, double, double,
				ivec shift, std::complex<double>,
				const symmetry& S, int sn, void* cld) {
    
    CylindricalChunkloopData* chunkloop_data = static_cast<CylindricalChunkloopData*>(cld);

    // index vectors for start and end of chunk
    ivec isS = S.transform(is, sn) + shift;
    ivec ieS = S.transform(ie, sn) + shift;

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

    // Build and record the chunk metadata
    ZRVector<std::size_t> chunk_shape{shape[0], shape[1]};
    ZRVector<std::size_t> downsampled_chunk_shape(to_downsampled_shape(chunk_shape, chunkloop_data -> downsampling_on_disk));
    CylindricalChunkloopData::sim_chunk_meta_t cur_meta(chunk_shape, downsampled_chunk_shape);    
    chunkloop_data -> sim_chunk_meta.emplace(std::make_pair(ichunk, cur_meta));

    // Setup field statistics tracker for this chunk
    chunkloop_data -> fstats.AddRegion(ichunk, chunk_shape);
  }

  // This is called for every simulation timestep and performs the saving of the Green's function.
  void eisvogel_saving_chunkloop(fields_chunk* fc, int ichunk, component cgrid, ivec is, ivec ie,
				 vec, vec, vec, vec, double, double,
				 ivec shift, std::complex<double> shift_phase,
				 const symmetry& S, int sn, void* cld) {
    
    CylindricalChunkloopData* chunkloop_data = static_cast<CylindricalChunkloopData*>(cld);    
    
    // Number of time slices before a new chunk is started
    std::size_t requested_chunk_size_t = 200;

    // Make sure the buffers are of the correct size for this simulation chunk
    ZRVector<std::size_t> spatial_chunk_shape(chunkloop_data -> sim_chunk_meta.at(ichunk).chunk_shape);
    TZRVector<std::size_t> field_slice_shape(1, spatial_chunk_shape);

    chunkloop_data -> field_absval_buffer.resize(spatial_chunk_shape);
    chunkloop_data -> field_buffer.resize(field_slice_shape);

    assert(is.z() >= 0);
    assert(is.r() >= 0);

    // Start index of this simulation chunk
    ZRIndexVector spatial_chunk_start_ind{
      (std::size_t)((is.z() - 1) / 2),
      (std::size_t)((is.r() - 1) / 2)
    };
    TZRIndexVector chunk_start_ind(chunkloop_data -> ind_time, spatial_chunk_start_ind);
      
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
      
      // fetch field components at child point ...
      data.update_values(idx);
      double E_z_val = data.values[0].real();
      double E_r_val = data.values[1].real();      
      double E_abs_val = std::sqrt(E_z_val * E_z_val + E_r_val * E_r_val);

      // ... and store them
      CylindricalChunkloopData::view_t field_elem = chunkloop_data -> field_buffer[cur_ind - chunk_start_ind];
      field_elem[0] = (scalar_t)E_r_val;
      field_elem[1] = (scalar_t)E_z_val;
      chunkloop_data -> field_absval_buffer[cur_spatial_ind - spatial_chunk_start_ind] = (scalar_t)E_abs_val;      
    }

    // Update statistics tracker
    chunkloop_data -> fstats.UpdateStatisticsForRegion(ichunk, chunkloop_data -> field_absval_buffer);

    // Truncate small field values to keep a certain maximum dynamic range
    limit_field_dynamic_range(ichunk, chunkloop_data -> field_buffer, chunkloop_data -> field_absval_buffer, chunkloop_data -> fstats,
			      chunkloop_data -> abs_min_field, chunkloop_data -> dynamic_range);

    // If requested, downsample the field buffer before it is registered in the output storage ...
    downsample_field(chunkloop_data -> field_buffer, chunkloop_data -> field_buffer_processed, chunkloop_data -> downsampling_on_disk);

    // ... and also compute the start index of the processed chunk (the index along the time axis remains untouched, any processing happens only in the spatial directions)
    ZRIndexVector spatial_processed_chunk_start_ind(to_downsampled_ind(spatial_chunk_start_ind, chunkloop_data -> downsampling_on_disk));
    TZRIndexVector processed_chunk_start_ind(chunkloop_data -> ind_time, spatial_processed_chunk_start_ind);
    
    using shape_t = CylindricalChunkloopData::chunk_t::shape_t;
    shape_t requested_storage_chunk_size(400);

    // Register the processed buffer
    if(chunkloop_data -> ind_time % requested_chunk_size_t == 0) {

      // Register this slice as the beginning of a new chunk ...
      TZRIndexVector start(0);
      auto simulation_chunk_storer = [&](const TZRIndexVector& storage_chunk_start, const TZRIndexVector& storage_chunk_end) {
	chunkloop_data -> field_chunk_buffer.resize(storage_chunk_end - storage_chunk_start);
	chunkloop_data -> field_chunk_buffer.fill_from(chunkloop_data -> field_buffer_processed,
						       storage_chunk_start,
						       storage_chunk_end,
						       start);
	chunkloop_data -> fstor.RegisterChunk(chunkloop_data -> field_chunk_buffer, processed_chunk_start_ind + storage_chunk_start);
      };
      IteratorUtils::index_loop_over_chunks(start, chunkloop_data -> field_buffer_processed.GetShape(), requested_storage_chunk_size, simulation_chunk_storer);
      
      // std::cout << "registering new chunk at start_ind = " << chunk_start_ind << " with shape = " << chunkloop_data -> field_buffer.GetShape() << std::endl;      
    }
    else {

      // ... or append it to an already-existing chunk along the outermost (time) direction
      TZRIndexVector start(0);
      auto simulation_chunk_appender = [&](const TZRIndexVector& storage_chunk_start, const TZRIndexVector& storage_chunk_end) {
	chunkloop_data -> field_chunk_buffer.resize(storage_chunk_end - storage_chunk_start);
	chunkloop_data -> field_chunk_buffer.fill_from(chunkloop_data -> field_buffer_processed,
						       storage_chunk_start,
						       storage_chunk_end,
						       start);
	chunkloop_data -> fstor.AppendSlice<0>(processed_chunk_start_ind + storage_chunk_start, chunkloop_data -> field_chunk_buffer);
      };
      IteratorUtils::index_loop_over_chunks(start, chunkloop_data -> field_buffer_processed.GetShape(), requested_storage_chunk_size, simulation_chunk_appender);
      
      // std::cout << "appending slice at start_ind = " << chunk_start_ind << " with shape = " << chunkloop_data -> field_buffer.GetShape() << std::endl;     
    }
  }  
} // end namespace meep

namespace GreensFunctionCalculator::MEEP {

  void CylindricalGreensFunctionCalculator::calculate_mpi_chunk(std::filesystem::path outdir, std::filesystem::path local_scratchdir,
								double courant_factor, double resolution, double pml_width, std::size_t downsampling_on_disk) {
    
    std::shared_ptr<meep::grid_volume> meep_gv = std::make_shared<meep::grid_volume>(meep::volcyl(m_geom.GetRMax(), m_geom.GetZMax() - m_geom.GetZMin(), resolution));
    std::shared_ptr<meep::structure> meep_s = std::make_shared<meep::structure>(*meep_gv, m_geom, meep::pml(pml_width), meep::identity(), 0, courant_factor);
    std::shared_ptr<meep::fields> meep_f = std::make_shared<meep::fields>(meep_s.get());
    
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
    
    scalar_t dynamic_range = 50;     // max. dynamic range of field strength to keep in the stored output
    scalar_t abs_min_field = 1e-20;  // absolute minimum of field to retain in the stored output
    
    // Prepare data container to pass to all MEEP callbacks
    TZRVector<std::size_t> init_field_buffer_shape(1);
    TZRVector<std::size_t> init_field_buffer_downsampled_shape(1);
    TZRVector<std::size_t> init_field_chunk_buffer_shape(1);
    ZRVector<std::size_t> init_field_absval_buffer_shape(1);
    CylindricalChunkloopData cld(0, darr, fstats, dynamic_range, abs_min_field, downsampling_on_disk,
				 init_field_buffer_shape, init_field_buffer_downsampled_shape,
				 init_field_absval_buffer_shape, init_field_chunk_buffer_shape);
    
    // Setup simulation run 
    meep_f -> loop_in_chunks(meep::eisvogel_setup_chunkloop, static_cast<void*>(&cld), meep_f -> total_volume());
    
    // Main simulation loop
    std::size_t stepcnt = 0;
    for(double cur_t = 0.0; cur_t <= m_t_end; cur_t += 0.1) {
      
      // Time-step the fields
      while (meep_f -> time() < cur_t) {
	meep_f -> step();
      }
      
      if(meep::am_master()) {
	std::cout << "Simulation time: " << meep_f -> time() << std::endl;
      }
      
      // Store the Green's function at the present timestep
      cld.ind_time = stepcnt++;
      meep_f -> loop_in_chunks(meep::eisvogel_saving_chunkloop, static_cast<void*>(&cld), meep_f -> total_volume());
    }
    
    // Swap axes 0 and 2 to go from TZR to RZT
    darr.SwapAxes<0, 2>();
    
    // Move Green's function to the requested output location
    darr.Move(outdir);  
  }
  
  CylindricalGreensFunctionCalculator::CylindricalGreensFunctionCalculator(CylinderGeometry& geom, const Antenna& antenna, scalar_t t_end) :
    m_t_end(t_end), m_geom(geom), m_antenna(antenna) {
    
    m_start_coords = std::make_shared<RZTCoordVector>((scalar_t)0.0, geom.GetZMin(), (scalar_t)0.0);
    m_end_coords = std::make_shared<RZTCoordVector>(geom.GetRMax(), geom.GetZMax(), t_end);
    
    std::cout << "Constructed geom with start = " << *m_start_coords << " and end = " << *m_end_coords << std::endl;
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
      darr.Import(indir);
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
							int cur_mpi_id, int number_mpi_jobs, const RZTVector<std::size_t>& requested_chunk_size, std::size_t overlap) {
    
    using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
    
    std::size_t cache_depth = 5;
    RZTVector<std::size_t> init_cache_el_shape(1);
    RZTVector<std::size_t> streamer_chunk_shape(stor::INFTY); streamer_chunk_shape[0] = 1; // serialize one radial slice at a time
    darr_t darr(indir, cache_depth, init_cache_el_shape, streamer_chunk_shape);
    
    std::cout << "loading array from " << indir << std::endl;
    std::cout << "now rechunking; have input array with shape " << darr.GetShape() << std::endl;
    std::cout << "this is rechunking job " << cur_mpi_id << " / " << number_mpi_jobs << " jobs" << std::endl;
    std::cout << "rechunking onto chunk_size = " << requested_chunk_size << " and overlap = " << overlap << std::endl;
    
    // Determine nonoverlapping rectangular partitions that are handled by each job; try to make them of similar size for optimal work sharing
    RZTVector<std::size_t> shape = darr.GetShape();
    RZTVector<std::size_t> partition_size = get_partition_size(shape, requested_chunk_size, number_mpi_jobs);
    
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
	darr.RebuildChunksPartial(partition_start, partition_end, requested_chunk_size, cur_rechunking_dir, overlap, SpatialSymmetry::Cylindrical<scalar_t>::boundary_evaluator);
	std::cout << "job " << cur_mpi_id << " done" << std::endl;
      }
      
      num_rechunking++;
    };
    RZTIndexVector start(0);
    IteratorUtils::index_loop_over_chunks(start, shape, partition_size, rechunker);
    
    // After all jobs are finished with rechunking their respective regions, merge all outputs into `outdir`
    meep::all_wait(); 
    
    if(meep::am_master()) {
      merge_mpi_chunks(outdir, rechunking_dirs);
    }
  }
  
  void CylindricalGreensFunctionCalculator::Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
						      double courant_factor, double resolution, double pml_width, std::size_t downsampling_on_disk) {
    
    int number_mpi_jobs = meep::count_processors();
    int job_mpi_rank = meep::my_rank();
    
    // Prepare an empty working directory in global scratch area
    std::filesystem::path global_workdir = global_scratchdir / "eisvogel";
    if(meep::am_master()) {
      if(std::filesystem::exists(global_workdir)) {
	std::filesystem::remove_all(global_workdir);
      }
      std::filesystem::create_directory(global_workdir);
    }
    
    // Output directory
    if(meep::am_master()) {
      if(std::filesystem::exists(outdir)) {
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
    calculate_mpi_chunk(job_outdir, local_scratchdir, courant_factor, resolution, pml_width, downsampling_on_disk);
    
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
    RZTVector<std::size_t> requested_chunk_size(400);
    std::size_t overlap = 2;
    rechunk_mpi(outdir, mergedir, global_workdir, job_mpi_rank, number_mpi_jobs, requested_chunk_size, overlap);
    
    meep::all_wait();
    
    // Clean up our temporary files
    if(meep::am_master()) {
      std::filesystem::remove_all(global_workdir);
      
      // Fourth stage: create the actual Green's function object
      using darr_t = typename SpatialSymmetry::Cylindrical<scalar_t>::darr_t;
      darr_t darr(outdir);
      RZTVector<std::size_t> darr_shape = darr.GetShape();
      RZTCoordVector step_size = (*m_end_coords - *m_start_coords) / (darr_shape.template as_type<scalar_t>() - 1);
      CylindricalGreensFunction(*m_start_coords, *m_end_coords, step_size, std::move(darr));    
    }
  }
}
