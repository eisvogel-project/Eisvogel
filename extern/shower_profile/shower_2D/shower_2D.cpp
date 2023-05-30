#include "shower_2D.h"
#include <iostream>
#include <array>
#include "Eisvogel/CoordUtils.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "constants.h"
#include "units.h"

using namespace CoordUtils;

showers::Shower2D::Shower2D(
                std::array<scalar_t, 3> pos,
                scalar_t en,
                scalar_t ze,
                scalar_t az,
                ChargeExcessProfile ce,
                scalar_t ce_scaling,
                environment::IceProfile &ice

){
    starting_position = pos;
    energy = en;
    zenith = ze;
    azimuth = az;
    charge_excess_profile = ce;
    charge_excess_profile_scaling = ce_scaling;
    ice_profile = ice;
}