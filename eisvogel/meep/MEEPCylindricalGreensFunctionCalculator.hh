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

    // dynamic_range ... max. dynamic range of field strength to keep in the stored output
    // abs_min_field ... absolute minimum of field to retain in the stored output
    void Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
		   double courant_factor = 0.5, double resolution = 24, double timestep = 0.1, double pml_width = 2.0, std::size_t downsampling_on_disk = 2,
		   scalar_t dynamic_range = 50, scalar_t abs_min_field = 1e-20, std::size_t chunk_overlap = 2, std::size_t chunk_size_linear = 400,
		   std::size_t rechunk_cache_depth = 12);
    
  private:

    // Main steps of the Green's function calculation
    void calculate_mpi_chunk(std::filesystem::path outdir, std::filesystem::path local_scratchdir, double courant_factor, double resolution, double timestep, double pml_width,
			     std::size_t downsampling_on_disk, scalar_t dynamic_range, scalar_t abs_min_field);
    static void merge_mpi_chunks(std::filesystem::path outdir, const std::vector<std::filesystem::path>& indirs);
    static void rechunk_mpi(std::filesystem::path outdir, std::filesystem::path indir, std::filesystem::path global_scratchdir, int cur_mpi_id, int number_mpi_jobs,
			    const RZTVector<std::size_t>& requested_chunk_size, std::size_t overlap, std::size_t cache_depth);
    
  private:
    
    using meep_chunk_ind_t = int;
    
    scalar_t m_t_end;
    CylinderGeometry& m_geom;
    Antenna& m_antenna;
    
    std::shared_ptr<RZTCoordVector> m_start_coords;
    std::shared_ptr<RZTCoordVector> m_end_coords;
  };
}
  
#include "MEEPCylindricalGreensFunctionCalculator.hxx"
