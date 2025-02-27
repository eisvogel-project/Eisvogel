#pragma once

#include <memory>
#include <filesystem>
#include <meep.hpp>
#include "Geometry.hh"
#include "Antenna.hh"

struct RZTCoordVector;

template <typename T>
struct RZTVector;

namespace GreensFunctionCalculator::MEEP {

  class CylindricalGreensFunctionCalculator {
    
  public:
    
    CylindricalGreensFunctionCalculator(CylinderGeometry& geom, Antenna& antenna, scalar_t t_end);

    CylindricalGreensFunctionCalculator(CylinderGeometry& geom, const CylinderRegion& request_to_store, Antenna& antenna, scalar_t t_end);

    // dynamic_range ... max. dynamic range of field strength to keep in the stored output
    // abs_min_field ... absolute minimum of field to retain in the stored output
    void Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
		   double courant_factor = 0.5, double resolution = 28, double timestep = 0.1, double pml_width = 1.0, std::size_t downsampling_on_disk = 2,
		   scalar_t dynamic_range = 100, scalar_t abs_min_field = 1e-20, std::size_t chunk_overlap = 2, std::size_t chunk_size_linear = 800,
		   std::size_t rechunk_cache_depth = 5);
    
  private:

    struct MPIChunkCalculationResult;
    
    // Main steps of the Green's function calculation
    MPIChunkCalculationResult calculate_mpi_chunk(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::size_t padding_requested,
						  double courant_factor, double resolution, double timestep, double pml_width,
						  std::size_t downsampling_on_disk, scalar_t dynamic_range, scalar_t abs_min_field);
    static void merge_mpi_chunks(std::filesystem::path outdir, const std::vector<std::filesystem::path>& indirs);
    static void rechunk_mpi(std::filesystem::path outdir, std::filesystem::path indir, std::filesystem::path global_scratchdir, int cur_mpi_id, int number_mpi_jobs,
			    const RZTVector<std::size_t>& requested_chunk_size, std::size_t overlap, std::size_t cache_depth);
    
  private:
    
    using meep_chunk_ind_t = int;
    
    scalar_t m_t_end;
    CylinderGeometry& m_geom;
    CylinderRegion m_request_to_store; // Spatial region that is requested to be stored
    Antenna& m_antenna;    
  };
}
  
#include "MEEPCylindricalGreensFunctionCalculator.hxx"
