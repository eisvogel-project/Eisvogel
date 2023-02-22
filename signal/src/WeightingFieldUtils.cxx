#include "WeightingFieldUtils.hh"
#include "Common.hh"
#include "IteratorUtils.hh"
#include "CoordUtils.hh"
#include "NDArray.hh"

namespace C = CoordUtils;

namespace WeightingFieldUtils {

  WeightingField<> CreateElectricDipoleWeightingField() {   
    
    // These will be arguments eventually
    CoordVector start_coords = C::MakeCoordVectorTRZ(0.0, 0.0, 0.0);
    CoordVector end_coords = C::MakeCoordVectorTRZ(1.0, 1.0, 1.0);
    DeltaVector delta = C::MakeCoordVectorTRZ(0.1, 0.1, 0.1); // step size
    
    // ==============

    CoordVector number_pts = (end_coords - start_coords) / delta;
    std::size_t pts_t = C::getT(number_pts), pts_r = C::getR(number_pts), pts_z = C::getZ(number_pts);

    auto E_r = [&](scalar_t t, scalar_t r, scalar_t z) -> scalar_t {
      return 0.0;
    };
    
    auto E_z = [&](scalar_t t, scalar_t r, scalar_t z) -> scalar_t {
      return 1.0;
    };
    
    auto E_phi = [&](scalar_t t, scalar_t r, scalar_t z) -> scalar_t {
      return 2.0;
    };
    
    // TODO: find a better way to make sure the ordering t, z, r etc is not messed up
    ScalarField3D<scalar_t> E_r_sampled({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> E_z_sampled({pts_t, pts_z, pts_r}, 0.0);
    ScalarField3D<scalar_t> E_phi_sampled({pts_t, pts_z, pts_r}, 0.0);
    
    IndexVector start_inds({0, 0, 0});
    IndexVector end_inds({pts_t, pts_z, pts_r});

    for(IndexCounter cnt(start_inds, end_inds); cnt.running(); ++cnt) {

      IndexVector ind = cnt.index();

      scalar_t t = C::getT(start_coords) + C::getTInd(ind) * C::getT(delta);
      scalar_t r = C::getR(start_coords) + C::getRInd(ind) * C::getR(delta);
      scalar_t z = C::getZ(start_coords) + C::getZInd(ind) * C::getZ(delta);
      
      E_r_sampled(ind) = E_r(t, r, z);
      E_z_sampled(ind) = E_z(t, r, z);
      E_phi_sampled(ind) = E_phi(t, r, z);
    }
    
    return WeightingField(std::move(E_r_sampled), std::move(E_z_sampled), std::move(E_phi_sampled),
			       std::move(start_coords), std::move(end_coords));
  }
}
