#pragma once

#include <fstream>
#include <filesystem>
#include "Eisvogel/Common.hh"
#include "Vector.hh"
#include "DistributedNDVecArray.hh"
#include "Interpolation.hh"

struct CylindricalGreensFunctionMetadata {

  CylindricalGreensFunctionMetadata();
  CylindricalGreensFunctionMetadata(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval);
  
  RZTCoordVector start_pos;
  RZTCoordVector end_pos;
  RZTCoordVector sample_interval;  
};

// Green's function for a cylindrically-symmetric geometry. It abstracts away the discreteness of the underlying data and presents itself
// as a continuous function that can be evaluated at an arbitrary space-time point.
template <class KernelT>
class CylindricalGreensFunction {
  
public:

  // To load an existing cylindrically-symmetric Green's function
  CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size = 4);

  // To create a Green's function from existing sampled data
  CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
			    Distributed_RZT_ErEz_Array&& data);

  ~CylindricalGreensFunction();
  
  // Calculate inner product with source current over the time interval [t_start, t_end) and accumulate into `result`
  void accumulate_inner(const XYZCoordVector& coords, scalar_t t_start, scalar_t t_end, scalar_t t_samp,
			const XYZTFieldVector& current, std::vector<scalar_t>::iterator result);
    
private:
   
  static std::filesystem::path move_path_from(Distributed_RZT_ErEz_Array&& data);

  void load_metadata();
  void save_metadata();
  
private:

  static constexpr std::size_t dims = 3;        // A cylindrically-symmetric Green's function is indexed as (R, Z, T) for m = 0 ...
  static constexpr std::size_t vec_dims = 2;    // ... and stores (Er, Ez) at each location
  static constexpr std::string_view m_meta_filename = "green.meta";     // Filename for Green's function metadata

  std::filesystem::path m_meta_path;
  CylindricalGreensFunctionMetadata m_meta;
  ChunkLibrary<NDVecArray, scalar_t, dims, vec_dims> m_lib;
};

#include "GreensFunction.hxx"
