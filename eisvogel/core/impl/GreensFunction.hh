#pragma once

#include <fstream>
#include <filesystem>
#include "Vector.hh"
#include "DistributedNDVecArray.hh"
#include "Quadrature.hh"
#include "NDVecArray.hh"
#include "Common.hh"
#include "Current.hh"

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

    using view_t = typename NDVecArray<scalar_t, dims, vec_dims>::view_t;
    using lib_t = ChunkLibrary<NDVecArray, scalar_t, dims, vec_dims>;
    using chunk_t = typename ChunkLibrary<NDVecArray, scalar_t, dims, vec_dims>::chunk_t;
  };

  enum OutOfBoundsBehavior {
    RaiseError, Ignore
  };  
}

// Green's function for a cylindrically-symmetric geometry. It abstracts away the discreteness of the underlying data and presents itself
// as a continuous function that can be evaluated at an arbitrary space-time point.
class CylindricalGreensFunction : private Green::DimTraits<CylindricalGreensFunction>::lib_t {

private:

  using lib_t = Green::DimTraits<CylindricalGreensFunction>::lib_t;
  using view_t = Green::DimTraits<CylindricalGreensFunction>::view_t;

public:

  using chunk_t = Green::DimTraits<CylindricalGreensFunction>::chunk_t;
  
public:

  // To load an existing cylindrically-symmetric Green's function
  CylindricalGreensFunction(std::filesystem::path path, std::size_t cache_size = 4);

  // To create a Green's function from existing sampled data
  CylindricalGreensFunction(const RZTCoordVector& start_pos, const RZTCoordVector& end_pos, const RZTCoordVector& sample_interval,
			    Distributed_RZT_ErEz_Array&& data);

  // Apply this Green's function to the line current segment `curr_seg` and accumulate the result in `signal`.
  // The signal is calculated starting from time `t_sig_start` with `num_samples` samples using `t_sig_samp` as sampling interval.
  template <class KernelT, typename ResultT = float, class QuadratureT = Quadrature::TrapezoidalRule>
  void apply_accumulate(const LineCurrentSegment& seg, scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples,
			std::vector<ResultT>& signal, 
			Green::OutOfBoundsBehavior oob_mode = Green::OutOfBoundsBehavior::RaiseError,
			scalar_t weight = 1.0, scalar_t max_itgr_step = 1.0);

  template <class KernelT>
  void fill_array(const RZTCoordVector& start_coords, const RZTCoordVector& end_coords, const RZTVector<std::size_t>& num_samples, chunk_t& array);

  RZTCoordVector start_coords();
  RZTCoordVector end_coords();

private:

  // Various specialized versions of `apply_accumulate`
  // General version that works for arbitrary line segments
  template <class KernelT, typename ResultT, class QuadratureT>
  void apply_accumulate_general(const LineCurrentSegment& seg, scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples,
				std::vector<ResultT>& signal, 
				Green::OutOfBoundsBehavior oob_mode = Green::OutOfBoundsBehavior::RaiseError,
				scalar_t weight = 1.0, scalar_t max_itgr_step = 1.0);

  // Optimized for short tracks that are shorter than `itgr_step`
  template <class KernelT, typename ResultT, class QuadratureT>
  void apply_accumulate_short_segment(const LineCurrentSegment& seg, scalar_t t_sig_start, scalar_t t_sig_samp, std::size_t num_samples,
				      std::vector<ResultT>& signal, 
				      Green::OutOfBoundsBehavior oob_mode = Green::OutOfBoundsBehavior::RaiseError,
				      scalar_t weight = 1.0);
  
public: // TODO: make this private, called by unit test at the moment
  
  // Calculate inner product with source current over the time interval [t_start, t_end) and accumulate into `result`
  template <class KernelT, typename ResultT>
  void accumulate_inner_product(const RZCoordVectorView coords, scalar_t t_start, scalar_t t_samp, std::size_t num_samples,
				const RZFieldVectorView source, std::vector<ResultT>::iterator result, scalar_t weight = 1.0f,
				Green::OutOfBoundsBehavior oob_mode = Green::OutOfBoundsBehavior::RaiseError);
    
private:

  // Converts cartesian (x, y, z) coordinates into cylindrical (r, z) coordinates
  void coord_cart_to_cyl(const XYZCoordVector& coords_cart, RZCoordVectorView coords_cyl);

  // Converts a field vector (A_x, A_y, A_z) defined at (x, y, z) into cylindrical field components (A_r, A_z)
  void field_cart_to_cyl(const XYZFieldVector& field_cart, const XYZCoordVector& coords_cart, RZFieldVectorView field_cyl);
  
  // Performs the inner product in cylindrical coordinates
  scalar_t inner_product(const RZFieldVectorView field, const RZFieldVectorView source);
  
  // Useful for converting from coordinates to storage indices
  RZCoordVector coords_to_index(const RZCoordVectorView rz_coords);
  RZTCoordVector coords_to_index(const RZTCoordVector& rzt_coords);
  
  // Used during (de)construction
  static std::filesystem::path move_path_from(Distributed_RZT_ErEz_Array&& data);
  void load_metadata();
  void save_metadata();
  
private:

  static constexpr std::size_t vec_dims = Green::DimTraits<CylindricalGreensFunction>::vec_dims;
  static constexpr std::string_view m_meta_filename = "green.meta";     // Filename for Green's function metadata

  CylindricalGreensFunctionMetadata m_meta;
  std::filesystem::path m_meta_path;
};

#include "GreensFunction.hxx"
