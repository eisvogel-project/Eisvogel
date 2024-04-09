#pragma once

#include <fstream>
#include <filesystem>
#include "Eisvogel/Common.hh"
#include "Vector.hh"
#include "DistributedNDVecArray.hh"

struct CylindricalGreensFunctionMetadata {

  // Can add more information here later, e.g. about the geometry that went into the Green's function
  
  CylindricalGreensFunctionMetadata();
  CylindricalGreensFunctionMetadata(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval);
  
  RZTCoordVector start_pos_rzt;
  RZTCoordVector end_pos_rzt;
  RZTCoordVector sample_interval_rzt;

  RZCoordVector start_pos_rz;
  RZCoordVector end_pos_rz;
  RZCoordVector sample_interval_rz;  
};

class CylindricalGreensFunction;

namespace Green {

  template <typename T>
  struct DimTraits;
  
  template <>
  struct DimTraits<CylindricalGreensFunction> {
    static constexpr std::size_t dims = 3;        // A cylindrically-symmetric Green's function is indexed as (R, Z, T) for m = 0 ...
    static constexpr std::size_t vec_dims = 2;    // ... and stores (Er, Ez) at each location

    using lib_t = ChunkLibrary<NDVecArray, scalar_t, dims, vec_dims>;
  };
}

// Green's function for a cylindrically-symmetric geometry. It abstracts away the discreteness of the underlying data and presents itself
// as a continuous function that can be evaluated at an arbitrary space-time point.
class CylindricalGreensFunction : private Green::DimTraits<CylindricalGreensFunction>::lib_t {

private:

  using lib_t = Green::DimTraits<CylindricalGreensFunction>::lib_t;
  
public:

  // To load an existing cylindrically-symmetric Green's function
  CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size = 4);

  // To create a Green's function from existing sampled data
  CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
			    Distributed_RZT_ErEz_Array&& data);
  
  // Calculate inner product with source current over the time interval [t_start, t_end) and accumulate into `result`
  template <class KernelT>
  void accumulate_inner_product(const RZCoordVector& coords, scalar_t t_start, scalar_t t_end, scalar_t t_samp,
				const XYZFieldVector& source, std::vector<scalar_t>::iterator result);
    
private:

  // Useful for converting from coordinates to storage indices
  RZCoordVector coords_to_index(const RZCoordVector& rz_coords);
  RZTCoordVector coords_to_index(const RZTCoordVector& rzt_coords);
  
  // Used during (de)construction
  static std::filesystem::path move_path_from(Distributed_RZT_ErEz_Array&& data);
  void load_metadata();
  void save_metadata();
  
private:
  
  static constexpr std::string_view m_meta_filename = "green.meta";     // Filename for Green's function metadata

  CylindricalGreensFunctionMetadata m_meta;
  std::filesystem::path m_meta_path;
};

#include "GreensFunction.hxx"
