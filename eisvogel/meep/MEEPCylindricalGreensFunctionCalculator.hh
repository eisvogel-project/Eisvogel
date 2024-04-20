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
    
    CylindricalGreensFunctionCalculator(CylinderGeometry& geom, const Antenna& antenna, scalar_t t_end);
    
    void Calculate(std::filesystem::path outdir, std::filesystem::path local_scratchdir, std::filesystem::path global_scratchdir,
		   double courant_factor = 0.5, double resolution = 12, double pml_width = 1.0);
    
  private:
    
    void calculate_mpi_chunk(std::filesystem::path outdir, std::filesystem::path local_scratchdir, double courant_factor, double resolution, double pml_width);
    static void merge_mpi_chunks(std::filesystem::path outdir, const std::vector<std::filesystem::path>& indirs);
    static void rechunk_mpi(std::filesystem::path outdir, std::filesystem::path indir, std::filesystem::path global_scratchdir, int cur_mpi_id, int number_mpi_jobs,
			    const RZTVector<std::size_t>& requested_chunk_size, std::size_t overlap);
    
  private:
    
    using meep_chunk_ind_t = int;
    
    scalar_t m_t_end;
    CylinderGeometry& m_geom;
    const Antenna& m_antenna;
    
    std::shared_ptr<RZTCoordVector> m_start_coords;
    std::shared_ptr<RZTCoordVector> m_end_coords;
  };
}
  
#include "MEEPCylindricalGreensFunctionCalculator.hxx"
