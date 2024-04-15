#pragma once

#include <filesystem>
#include "Vector.hh"

namespace GreensFunctionUtils {

  void CreateElectricDipoleGreensFunction(std::filesystem::path gf_path, const RZCoordVector& start_coords, const RZCoordVector& end_coords, scalar_t t_end, scalar_t ior,
					  scalar_t filter_t_peak, unsigned int filter_order, scalar_t r_min,
					  scalar_t os_factor = 10, std::size_t max_pts_in_chunk = 400);

}

#include "GreensFunctionUtils.hxx"
